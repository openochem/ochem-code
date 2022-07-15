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

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.HibernateException;
import org.hibernate.criterion.Projections;

import qspr.Globals;
import qspr.entities.Property;
import qspr.entities.User;

@XmlRootElement
public class UserProfile {
	public User user;
	
	@XmlElementWrapper(name = "uploaded-properties")
	@XmlElement(name = "property")
	public List<Property> uploadedProperties;
	
	public Long pendingTasksToBeDeleted;
	public Long modelsToBeDeleted;
	public boolean firstLogin = false;
	
	public double bonusPoints;
	
	public static UserProfile create(User user) throws HibernateException, Exception {
		UserProfile profile = new UserProfile();
		profile.user = user;
		profile.uploadedProperties = new ArrayList<Property>();
		
//		Criteria c = Globals.session().createCriteria(ExperimentalProperty.class, "e")
//				.add(Restrictions.eq("introducer", user))
//				.setProjection(Projections.projectionList().add(Projections.groupProperty("e.property")).add(Projections.count("id"), "cnt"))
//				.addOrder(Order.desc("cnt"));
//		
//		if (!Globals.isSuperUser())
//			c.add(Restrictions.ge("rights", Globals.RIGHTS_WRITE));
//			
//		profile.uploadedProperties = new ArrayList<Property>();
//		List<Object[]> rows = c.list();
//		for (Object[] row : rows) 
//		{
//			profile.uploadedProperties.add((Property)row[0]);
//			((Property)row[0]).count = (Long) row[1];
//		}
		
		if (user.equals(Globals.myself()))
		{
			PendingTaskFilter pFilter = new PendingTaskFilter();
			pFilter.toBeDeleted = true;
			profile.pendingTasksToBeDeleted = (Long) pFilter.createCriteria().setProjection(Projections.count("id")).uniqueResult();
			
			ModelFilter filter = new ModelFilter();
			filter.toBeDeleted = true;
			profile.modelsToBeDeleted = (Long) new ModelOperation().getFilterCriteria(filter).setProjection(Projections.count("id")).uniqueResult();
		}
		
		return profile;
	}
}
