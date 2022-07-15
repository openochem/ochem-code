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

import com.eadmet.batchupload.main.RecordStub.NameValueStub;
import com.eadmet.batchupload.main.RecordStub.PropertyStub;

@XmlRootElement
public class PropertyRemapping extends AbstractPropertyRemapping
{
	public SchemaList<ConditionRemapping> conditions = new SchemaList<ConditionRemapping>().setGenericClass(ConditionRemapping.class);
	
	@Override
	public boolean valid(boolean warningIsInvalid)
	{
		if (!super.valid(warningIsInvalid))
			return false;
		for (ConditionRemapping c : conditions)
			if (!c.valid(warningIsInvalid))
				return false;
		return true;
	}
	
	public String toString()
	{
		return super.toString()+":[options: "+options+"; units: "+units+"; conditions: "+conditions+"]";
	}
	
	public void updateWith(PropertyStub stub)
	{
		super.updateWith(stub);
		
		for (NameValueStub conditionStub : stub.conditions) 
			conditions.getByName(conditionStub.name.value).updateWith(conditionStub);

	}
}