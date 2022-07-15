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

package qspr.exception;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.entities.Unit;

import com.eadmet.exceptions.UserFriendlyException;

@XmlRootElement(name = "unit-conversion-exception")
public class UnitConversionException extends UserFriendlyException
{
	private static final long serialVersionUID = 1L;

	@XmlElement
	public Unit unit;
	
	@XmlElement
	public Unit defaultUnit;
	
	@XmlTransient
	public boolean toDefault;
	
	public UnitConversionException()
	{
		
	}
	
	public String getMessage()
	{
		if (!toDefault)
			return "Cannot convert from default unit, " + defaultUnit + " to " + unit + " for categgory "+ unit.category;
		else
			return "Cannot convert from " + unit + " to default unit, " + defaultUnit + " for categgory "+ unit.category;
			
	}
	
	public UnitConversionException(Unit unit, boolean toDefault)
	{
		defaultUnit = unit.category.getDefaultUnit();
		this.toDefault = toDefault;
		this.unit = unit;
	}
}
