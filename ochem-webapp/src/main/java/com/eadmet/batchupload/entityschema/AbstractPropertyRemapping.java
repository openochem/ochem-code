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

package com.eadmet.batchupload.entityschema;

import javax.xml.bind.annotation.XmlRootElement;

import qspr.entities.PropertyValue;

import com.eadmet.batchupload.main.RecordStub.NameValueStub;

@XmlRootElement
public class AbstractPropertyRemapping extends RemappedValue
{
	public SchemaList<RemappedValue> options = new SchemaList<RemappedValue>().setGenericClass(RemappedValue.class);
	public SchemaList<RemappedValue> units = new SchemaList<RemappedValue>().setGenericClass(RemappedValue.class);
	
	@Override
	public boolean valid(boolean warningIsInvalid)
	{
		if (!super.valid(warningIsInvalid))
			return false;
		for (RemappedValue unit : units)
			if (!unit.valid(warningIsInvalid))
				return false;
		for (RemappedValue option: options)
			if (!option.valid(warningIsInvalid))
				return false;
		return true;
	}
	
	public void updateWith(NameValueStub stub)
	{
		super.updateWith(stub.name);
		
		Double d = null;
		try
		{
			Object[] vals = PropertyValue.parseTextIntoPredicateAndValues(stub.value.getFirst().value.value);
			d = (Double)vals[1];
		} catch (Exception e)
		{
			////
		}
		
		if (d != null) // Numeric value
		{
			if (stub.value.getFirst().unit.size() > 0)
			{
				RemappedValue rv = units.getByName(stub.value.getFirst().unit.getFirst().value);
				rv.updateWith(stub.value.getFirst().unit.getFirst());
				if (rv.minValue == null || d < rv.minValue)
					rv.minValue = d;
				if (rv.maxValue == null || d > rv.maxValue)
					rv.maxValue = d;
			}
			
		} else
			options.getByName(stub.value.getFirst().value.value).updateWith(stub.value.getFirst().value);
	}
	
	public String toString()
	{
		return super.toString()+":[options: "+options+"; units: "+units+"]";
	}
}