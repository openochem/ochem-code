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

import org.hibernate.Criteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.Session;
import qspr.entities.User;
import qspr.workflow.utils.QSPRConstants;

// Midnighter on Mar 22, 2012
public class BasketPeer 
{
	public static Criteria getListCriteria(BasketFilter filter)
	{
		Session session = Globals.userSession();
		Criteria basketCriteria = Globals.session().createCriteria(Basket.class, "basket");
		Disjunction disjunction = Restrictions.disjunction();

		if (filter.id != null)
			basketCriteria.add(Restrictions.eq("id", filter.id));

		if (session.user != null)
			if (session.user.group == null || !filter.showGroupBaskets)
				disjunction.add(Restrictions.eq("user", session.user));
			else
			{
				// Allow to see baskets of my group
				basketCriteria.createAlias("user", "u");
				disjunction.add(Restrictions.eq("u.group", session.user.group));
			}
		else
			disjunction.add(Restrictions.eq("session", session));

		if (filter.showPublicBaskets)
		{
			User published = Repository.user.getById(QSPRConstants.PUBLISHER_ID);
			if(filter.id == null)disjunction = Restrictions.disjunction(); // only public sets will be shown
			disjunction.add(Restrictions.eq("user", published));
		}

		if (filter.name != null)
			basketCriteria.add(Restrictions.like("name", "%" + filter.name + "%"));

		if (!filter.showSystemBaskets)
			basketCriteria.add(Restrictions.eq("basketType",0L));

		basketCriteria.add(disjunction);
		basketCriteria.addOrder(Order.desc("basketType"));
		basketCriteria.addOrder(Order.desc("id"));

		if (filter.mapping2 != null)
		{
			basketCriteria.createAlias("entries", "e");
			basketCriteria.createAlias("e.ep", "ep");
			basketCriteria.createAlias("ep.molecule", "mol");
			basketCriteria.add(Restrictions.eq("mol.mapping2", filter.mapping2));
		}

		return basketCriteria;
	}
}
