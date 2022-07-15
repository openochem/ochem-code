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

package qspr.util;

import java.util.List;

import qspr.entities.Basket;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Property;

public class BasicRecordMapper 
{
	private List<Property> properties;
	
	public Long getClass(ExperimentalProperty ep)
	{
		// TODO novserj: Check if it's correct way to do it
		if (ep.property.parent != null)
			return Long.valueOf(properties.indexOf(ep.property.parent));
		else
			return Long.valueOf(properties.indexOf(ep.property));
	}
	
	public Long getClass(Property property)
	{
		//FIXME: Check if it's correct way to do it
		if (property.parent != null)
			return Long.valueOf(properties.indexOf(property.parent));
		else
			return Long.valueOf(properties.indexOf(property));
	}
	
	public int getRowsSize()
	{
		return this.properties.size();
	}
	
	public BasicRecordMapper(Basket basket)
	{
		properties = basket.getProperty();
	}
}
