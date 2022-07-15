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

import java.util.ArrayList;
import java.util.List;

import org.hibernate.Criteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.User;

/**
 * Filter/query service for models
 * TODO: All model queries should be refactored to use this class
 * 
 * @author midnighter
 */
public class ModelsFilter
{
	public boolean pendingTasks = false;
	public List<Basket> trainingSets = new ArrayList<Basket>();
	public boolean featuredOnly = false;
	public boolean publicOnly = true;

	public Criteria createCriteria()
	{
		Criteria criteria = Globals.session().createCriteria(Model.class);
		criteria.createAlias("session", "sess");
		Disjunction authCriteria = Restrictions.disjunction();
		authCriteria.add(Restrictions.eq("session", Globals.userSession()));
		if (publicOnly)
			authCriteria.add(Restrictions.eq("published", new Boolean(true)));
		if (featuredOnly)
			authCriteria.add(Restrictions.isNotNull("featuredName"));
		if (Globals.userSession().user != null) {
			User user = Globals.userSession().user;
			//authCriteria.add(Restrictions.eq("sess.user", user)); // instead of this we have:
			if (user.group == null)
				authCriteria.add(Restrictions.eq("sess.user", user));
			else {
				criteria.createAlias("sess.user", "u");
				authCriteria.add(Restrictions.eq("u.group", user.group));
			}

		}
		criteria.add(authCriteria);

		if (!trainingSets.isEmpty())
		{
			if (trainingSets.size() == 1)
				criteria.add(Restrictions.eq("trainingSet", trainingSets.get(0)));
			else
				criteria.add(Restrictions.in("trainingSet", trainingSets));
		}

		if (pendingTasks)
			criteria.add(Restrictions.isNotNull("taskId"));
		else
			criteria.add(Restrictions.isNull("taskId"));

		return criteria;
	}
}
