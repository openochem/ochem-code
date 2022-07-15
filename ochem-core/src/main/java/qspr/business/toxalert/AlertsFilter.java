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

package qspr.business.toxalert;

import javax.servlet.http.HttpServletRequest;

import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.Property;
import qspr.entities.SubstructureAlert;
import qspr.entities.User;
import qspr.util.AccessChecker;
import qspr.util.RequestParser;

/**
 * Allows filtering structural alerts by various criteria
 * 
 * @author midnighter
 *
 */


public class AlertsFilter 
{
	public Long publicationId;
	public Long endpointId;
	public Long alertId;
	public String name;
	public boolean selectedOnly;
	public String order;
	public boolean includeFolders = false;
	public Long introducerId;
	
	public boolean approvedOnly = false;
	public boolean awaitingApproval = false;
	
	/**
	 * Add the filters to a criteria
	 * @param criteria - the hibernate criteria to modify
	 * @return the same criteria (for method chaining)
	 * @throws Exception
	 */
	public Criteria filterCriteria(Criteria criteria) throws Exception
	{
		if (criteria == null)
			criteria = Globals.session().createCriteria(SubstructureAlert.class);
		if (alertId != null)
			criteria.add(Restrictions.eq("id", alertId));
		if (endpointId != null)
			criteria.add(Restrictions.eq("property",
					Property.getById(endpointId)));
		if (publicationId != null)
			criteria.add(Restrictions.eq("article",
					Article.getById(publicationId)));
		if (name != null)
			criteria.add(Restrictions.like("name", "%"+name+"%"));
		
		if (introducerId != null)
			criteria.add(Restrictions.eq("introducer",
					User.getById(introducerId)));
		
		if (approvedOnly)
			criteria.add(Restrictions.eq("approved", true));
		
		if (awaitingApproval)
			criteria.add(Restrictions.or(Restrictions.eq("approved", false), Restrictions.isNull("approved")));
		
		if (!includeFolders)
			criteria.add(Restrictions.eq("folder", false));
		
		if (selectedOnly)
			if (!Globals.userSession().selectedAlerts.isEmpty())
				criteria.add(Restrictions.in("id", Globals.userSession().selectedAlerts));
			else
				criteria.add(Restrictions.eq("id", -1L)); // the selection is empty, show nothing
		
		User accessingUser = null;
		if (Globals.userSession() != null)
			accessingUser = Globals.userSession().user;
		
		AccessChecker.addAccessRestrictions(criteria, Globals.RIGHTS_FREELY_AVAILABLE, accessingUser, null, !approvedOnly);
		
		if (order != null)
			criteria.addOrder(Order.asc(order));
		
		return criteria;
	}
	
	/**
	 * Load filters from HTTP-based UI
	 */
	public void parseUI(HttpServletRequest request)
	{
		RequestParser parser = new RequestParser(request);
		if (parser.assertParam("id"))
			alertId = parser.getLongParam("id");
		if (parser.assertParam("alert-id"))
			alertId = parser.getLongParam("alert-id");
		if (parser.assertParam("property"))
			endpointId = parser.getLongParam("property");
		if (parser.assertParam("article"))
			publicationId = parser.getLongParam("article");
		if (parser.assertParam("alert-name"))
		{
			name = parser.getParam("alert-name");
			if (name.matches("TA[0-9]+")) // Special case - alert ID provided
			{
				alertId = Long.valueOf(name.substring(2));
				name = null;
			}
		}
		selectedOnly = parser.assertParam("selected-only");
		approvedOnly = parser.assertParam("approved-only");
		awaitingApproval = parser.assertParam("awaiting-approval");
	}
	
	/**
	 * Create a new criteria based on the filters
	 * @return the new criteria
	 * @throws Exception
	 */
	public Criteria filterCriteria() throws Exception
	{
		return filterCriteria(Globals.session().createCriteria(SubstructureAlert.class));
	}
	
	
}
