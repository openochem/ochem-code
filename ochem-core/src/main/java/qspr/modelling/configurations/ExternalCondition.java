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

package qspr.modelling.configurations;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.NumericalValueStandardizer;

import qspr.dao.Repository;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.Unit;
import qspr.metaserver.util.ShortCondition;

@XmlRootElement(name = "external-descriptor")
public class ExternalCondition extends ShortCondition
{
	private static final long serialVersionUID = 1L;

	private transient Property property;
	private transient Unit unit;

	public Property getProperty()
	{
		if (property != null)
			return property;
		return property = Repository.property.getPropertyById(id);
	}

	public Unit getUnit()
	{
		if (unit != null)
			return unit;
		return unit = Repository.unit.getUnitById(unitId);
	}

	public ExternalCondition()
	{
		// To make JAXB happy
	}

	public ExternalCondition(long id)
	{
		super(id);
	}

	public ExternalCondition(long id, String defaultOption)
	{
		super(id);
		this.defaultValue = Double.valueOf(getProperty().getOptionByName(defaultOption).id);
	}

	public ExternalCondition(ShortCondition desc) {
		super(desc);
	}


	public String toString()
	{
		return getProperty().getName();
	}

	private static String toString(Property property, Long modelId) {
		return " for property \"" + property.getName() + "\" for model: \"" + Repository.model.getByPublicId(modelId).name + "\"  with public id: " + modelId;
	}

	private static String getValueName(Property property, Double val) {
		PropertyOption o = PropertyOption.getById(val.longValue());
		if (property.isQualitative())return o != null? o.name : "NOT FOUND";
		return NumericalValueStandardizer.getSignificantDigits(val);
	}

	/**
	 *  Merges external oldOnes conditions (coming from DataSet or previous models) with newOnes (coming from additional models)
	 *  If there is no conflict, new ones will be added or merged
	 * @param oldOnes
	 * @param newOnes
	 * @param publicId
	 * @return
	 */

	public static List<ExternalCondition> checkAndMerge(List<ExternalCondition> oldOnes, List<ShortCondition> newOnes, Long modelId) {

		String add = "\n The selected(default) conditions as well as those of all models should be exactly the same.";

		for(ShortCondition a: newOnes) {

			Property propertyA = Repository.property.getPropertyById(a.id);
			Unit unitA = Repository.unit.getUnitById(a.unitId);

			boolean found = false;
			if(oldOnes != null)
				for(ExternalCondition b: oldOnes) {

					if(a.id == b.id) {
						found = true;
						if(!NumericalValueStandardizer.getSignificantDigits(a.defaultValue).equals(NumericalValueStandardizer.getSignificantDigits(b.defaultValue.doubleValue())))
							throw new UserFriendlyException("Default value: \"" + getValueName(propertyA, a.defaultValue) + "\" != \"" + 
									getValueName(b.getProperty(), b.defaultValue) + "\" " + toString(propertyA, modelId) + add);

						if(a.unitId != b.unitId)
							throw new UserFriendlyException("Units \"" + unitA  + "\" != \"" + b.getUnit() + "\"" + toString(propertyA, modelId) + add);

						if(a.optionsMergings != null && a.optionsMergings.size() > 0)
							throw new UserFriendlyException("Current implementation does not work with option merging.");

						if(b.optionsMergings != null &&  b.optionsMergings.size() > 0)
							throw new UserFriendlyException("Current implementation does not work with option merging.");
					}

				}

			if(!found) {
				if(oldOnes == null)oldOnes = new ArrayList<ExternalCondition>();
				oldOnes.add(new ExternalCondition(a));
			}
		}

		return oldOnes;
	}
}
