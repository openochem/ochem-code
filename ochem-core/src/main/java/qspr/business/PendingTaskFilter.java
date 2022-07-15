/* Copyright (C) 2022 BIGCHEM GmbH <info@bigchem.de>
 *
 * Contact: info@bigchem.de
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License (AGPL)
 * as published by the Free Software Foundation; either version 3.0
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the Affero GNU General Public License for more details.
 *
 * You should have received a copy of the Affero GNU Lesser General Public License
 * along with this program; If not, see <https://www.gnu.org/licenses/>. 
 */

package qspr.business;

import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.Criteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.PendingTask;
import qspr.entities.PendingTask.TaskType;

// Encapsulate the filter for pending tasks
// Midnighter on Jun 11, 2012

@XmlRootElement
public class PendingTaskFilter 
{
	public Long setID;
	public Long modelID;
	public String status;
	public TaskType taskType;
	public String taskName;

	public Long articleID;

	public boolean publishedOnly;
	public boolean ownTasksOnly = true;
	public boolean showGroupModels = false;

	public boolean toBeDeleted;

	public boolean unpublishedTasksOfOtherUsers = false;


	/**
	 * It was a quick fix to add visibility of Pending tasks to group
	 //TODO Refactor this code!
	 * @param c
	 * @param showGroupModels
	 */

	public static void addAuthRestrictions(Criteria c, boolean showGroupModels)
	{
		c.createAlias("session", "sess");
		Disjunction authCriteria = Restrictions.disjunction();
		authCriteria.add(Restrictions.eq("sess.id", Globals.userSession().id));
		if (Globals.userSession().user != null)
			if (Globals.userSession().user.group != null && showGroupModels)
			{
				c.createAlias("sess.user", "u", Criteria.LEFT_JOIN);
				authCriteria.add(Restrictions.eq("u.group", Globals.userSession().user.group));
			}
			else
				authCriteria.add(Restrictions.eq("sess.user", Globals.userSession().user));
		c.add(authCriteria);
	}


	public Criteria createCriteria() throws Exception
	{
		Criteria c = Globals.session().createCriteria(PendingTask.class);
		c.createAlias("model", "mdl", Criteria.LEFT_JOIN);

		if (ownTasksOnly && Globals.userSession().user.group != null && showGroupModels)
			addAuthRestrictions(c,showGroupModels);
		else{
			c.createAlias("session", "sess", Criteria.LEFT_JOIN);

			if (!unpublishedTasksOfOtherUsers || !Globals.isSuperUser())
			{
				if (!ownTasksOnly || !publishedOnly)
				{
					Disjunction disj = Restrictions.disjunction();
					disj.add(Restrictions.eq("published", true));
					if (!Globals.isGuestUser())
						disj.add(Restrictions.eq("sess.user", Globals.userSession().user));
					else
						disj.add(Restrictions.eq("session", Globals.userSession()));
					c.add(disj);
				}

				if (ownTasksOnly)
					if (!Globals.isGuestUser())
						c.add(Restrictions.eq("sess.user", Globals.userSession().user));
					else
						c.add(Restrictions.eq("session", Globals.userSession()));
			}
		}
		if (toBeDeleted)
			c.add(Restrictions.eqOrIsNull("mdl.published", false))
			.add(Restrictions.isNotNull("mdl.deleteAfterDays"));

		if (publishedOnly)
			c.add(Restrictions.eq("published", true));

		if (articleID != null)
			c.add(Restrictions.eq("article", Article.getById(articleID)));

		if (setID != null)
			c.createAlias("mdl.trainingSet", "ts").add(Restrictions.eq("ts.id", setID));

		if (modelID != null)
			c.add(Restrictions.eq("mdl.id", modelID));

		c.addOrder(Order.desc("timePosted"));

		if (status != null)
		{
			String[] statuses = status.split(",");
			if (statuses.length == 1)
				c.add(Restrictions.eq("status", status));
			else
				c.add(Restrictions.in("status", statuses));
		}

		if (taskType != null)
		{
			c.add(Restrictions.eq("type", taskType));
		}

		if (taskName != null)
		{
			if (taskName.matches("[0-9]+"))
				c.add(Restrictions.or(
						Restrictions.like("mdl.name", "%"+taskName+"%"), 
						Restrictions.eq("mdl.publicId", Long.valueOf(taskName))));
			else	
				c.add(Restrictions.or(
						Restrictions.like("mdl.name", "%"+taskName+"%"),
						Restrictions.like("name", "%"+taskName+"%"),
						Restrictions.like("setDescription", "%"+taskName+"%")));
		}

		return c;
	}
}