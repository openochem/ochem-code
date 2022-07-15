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

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hibernate.Criteria;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.business.ExperimentalRecordsFilter.ApprovalStatus;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Property;
import qspr.entities.PropertyOption;

import com.eadmet.exceptions.UserFriendlyException;

/**
 * A class for querying experimental records.
 * Its under construction.. 
 * A challenging project is to delayer the spaghetti code from ExperimentalPropertyBrowser
 * 
 * @author midnighter
 *
 */
public class ExperimentalRecordsPeer
{
	public static Criteria getCriteria(ExperimentalRecordsFilter filter)
	{
		Criteria c = Globals.session().createCriteria(ExperimentalProperty.class);

		if (filter.propertyId != null)
		{
			Property p = Property.getById(filter.propertyId);
			c.add(Restrictions.eq("property", p));
			if (filter.propertyOptionName != null)
				c.add(Restrictions.eq("option", p.getOption(filter.propertyOptionName)));
		}

		if (filter.basketId != null)
		{
			c.createAlias("basketEntries", "be");
			c.createAlias("be.basket", "b");
			c.add(Restrictions.eq("b.id", filter.basketId));
		}

		if (filter.basketId == null && filter.propertyId == null && filter.propertyOptionId == null)
			throw new UserFriendlyException("Too generic filter");

		if (filter.propertyOptionId != null)
			c.add(Restrictions.eq("option", PropertyOption.getById(filter.propertyOptionId)));


		ExperimentalProperty.addAccessRestrictions(c, Globals.RIGHTS_FREELY_AVAILABLE, Globals.userSession().user, filter.basketId != null, filter.dataFromOtherUsers != ApprovalStatus.APPROVED_ONLY);

		return c;
	}

	@SuppressWarnings("unchecked")
	public static Map<Long, Integer> recordsToMoleculeIds(List<Long> recordIds) {
		Map<Long, Integer> map = new HashMap<Long, Integer>();
		int processed = 0;
		while (processed < recordIds.size()) {
			int size = Math.min(1000, recordIds.size() - processed);
			List<Long> batchIds = recordIds.subList(processed, processed + size);
			processed += size;
			List<Object[]> rows = Globals.session().createCriteria(ExperimentalProperty.class)
					.createAlias("molecule", "mol")
					.createAlias("mol.mapping2", "mp2")
					.add(Restrictions.in("id", batchIds))
					.setProjection(Projections.projectionList().add(Projections.property("id")).add(Projections.property("mp2.id"))).list();
			for (Object[] row : rows) {
				map.put((Long)row[0], (Integer) row[1]);
			}
		}

		return map;
	}
}
