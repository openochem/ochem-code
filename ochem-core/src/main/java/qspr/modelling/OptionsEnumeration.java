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

package qspr.modelling;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import qspr.entities.Basket;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.util.BasicRecordMapper;

import com.eadmet.exceptions.UserFriendlyException;

public class OptionsEnumeration implements Serializable
{
	private static final long serialVersionUID = 1L;
	
	public Map<Long, Long> optionsMapping = new HashMap<Long, Long>();
	public List<Integer> numOfOptions = new ArrayList<Integer>();

	public static OptionsEnumeration enumeratePropertyOptions(Basket trainingSet, BasicRecordMapper mapper)
	{
		OptionsEnumeration en = new OptionsEnumeration();

		if (trainingSet == null)
			throw new UserFriendlyException("Empty training set provided");

		trainingSet.countRecordsByProperties(false);
		List<Property> properties = trainingSet.propertyUsed;

		Map<Long, String> optionNames = new HashMap<Long, String>();
		List<List<Long>> optionsByClass = new ArrayList<List<Long>>();

		for (int i = 0; i < mapper.getRowsSize(); i++)
			optionsByClass.add(new ArrayList<Long>());

		for(Property property : properties){

			if(!property.isQualitative())continue;
			
			List <PropertyOption> options = trainingSet.getPropertyOptions(property);

			for (PropertyOption option : options)
			{
				int classNum = mapper.getClass(property).intValue();

				if (!optionsByClass.get(classNum).contains(option.id))
					optionsByClass.get(classNum).add(option.id);

				optionNames.put(option.id, option.name);
			}
		}

		for (List<Long> classOptions : optionsByClass)
			Collections.sort(classOptions);

		for (List<Long> classOptions : optionsByClass)
			for (Long option : classOptions)
			{
				en.optionsMapping.put(option, Long.valueOf(classOptions.indexOf(option)));
				System.out.printf("<%s> = %d\n", optionNames.get(option), classOptions.indexOf(option));
			}

		en.numOfOptions = new ArrayList<Integer>();

		for (int i = 0; i < mapper.getRowsSize(); i++)
			en.numOfOptions.add(Math.max(1, optionsByClass.get(i).size()));
		return en;
	}
}
