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

package com.eadmet.useractions;

import javax.xml.bind.annotation.XmlRootElement;

import qspr.entities.Property;

@XmlRootElement
public class PropertyAction extends AbstractUserAction
{
	Property property;
	PropertyActionType type;
	
	@Override
	public String getLogLine()
	{
		switch (type)
		{
		case CREATE:
			return String.format(" has created a new property \"%s\"", property.getName());
		case DELETE:
			return String.format(" has deleted a property \"%s\"", property.getName());
		case RENAME:
			return String.format(" has created a property to \"%s\"", property.getName());
		}
		return null;
	}
	
	public PropertyAction()
	{
		
	}
	
	public PropertyAction(Property property, PropertyActionType type)
	{
		this.property = new Property();
		this.property.id = property.id;
		this.property.setName(property.getName());
		this.type = type;
	}
	
	public static enum PropertyActionType {CREATE, DELETE, RENAME} ; 
	
}
