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

package qspr.dao;

import java.util.List;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.ConditionSet;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.PropertyValue;
import qspr.interfaces.ProvidedConditions;
import qspr.metaserver.util.ShortCondition;
import qspr.modelling.configurations.ExternalCondition;

public class PropertyDAOImpl implements PropertyDAO {

	@Override
	@SuppressWarnings("unchecked")
	public Property getProperty(String name, boolean createIfMissing) {
		List<Property> props =	Globals.session().createCriteria(Property.class)
				.add(Restrictions.eq("shortName", Property.shortName(name)))
				.list();

		if (props.size() > 0)
			return props.get(0);

		if (!createIfMissing)
			return null;

		Property property = new Property();
		property.setName(name);
		return property;
	}

	@Override
	public Property getPropertyById(Long propertyId) 
	{
		return (Property)Globals.session().get(Property.class, propertyId);
	}

	/**
	 *  Conditions can be missing for some data points: 
	 *  We thus create a default set of them using configurations which supports them
	 * @param conf
	 * @return
	 */

	@Override
	public ConditionSet createDefaultConditions(ProvidedConditions conf, boolean startTransactions) 
	{
		if(startTransactions)Globals.startAllTransactions();

		List<ShortCondition> conditions = conf.getConditions();

		ConditionSet set = new ConditionSet();

		for(ShortCondition dd :conditions) {

			ExternalCondition d = dd instanceof ExternalCondition ? (ExternalCondition) dd : new ExternalCondition(dd);
			PropertyValue pv = new PropertyValue();
			pv.property = d.getProperty();
			if (pv.property.isNumeric()) {
				pv.value = d.defaultValue;
				pv.unit = d.getUnit();
			}
			else
				pv.option = (PropertyOption) Globals.session().get(PropertyOption.class, d.defaultValue.longValue());
			set.values.add(pv);
		}
		if(startTransactions)Globals.commitAllTransactions();

		return set;
	}

}
