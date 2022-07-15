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

import java.util.Map;

import javax.xml.bind.annotation.XmlRootElement;

import com.eadmet.batchupload.main.RecordStub;
import com.eadmet.batchupload.main.RecordStub.ColumnValue;

@XmlRootElement(name="remapping")
public class EntitiesRemapping
{
	public String fileSchemaMD5;
	public String md5;
	public Integer records = 0;
	
	public SchemaList<PropertyRemapping> properties = new SchemaList<PropertyRemapping>().setGenericClass(PropertyRemapping.class);
	public SchemaList<RemappedValue> articles = new SchemaList<RemappedValue>().setGenericClass(RemappedValue.class);
	public SchemaList<RemappedValue> baskets = new SchemaList<RemappedValue>().setGenericClass(RemappedValue.class);
	
	public boolean valid(boolean warningIsInvalid)
	{
		for (PropertyRemapping p : properties)
			if (!p.valid(warningIsInvalid))
				return false;
		for (RemappedValue article : articles)
			if (!article.valid(warningIsInvalid))
				return false;
		for (RemappedValue basket : baskets)
			if (!basket.valid(warningIsInvalid))
				return false;
		return true;
	}
	
	public void updateWith(RecordStub stub)
	{
		if (stub.property.name == null)
			return;

		PropertyRemapping property = properties.getByName(stub.property.name.value);
		property.updateWith(stub.property);
		
		if (stub.article.size() > 0)
			articles.getByName(stub.article.getFirst().value).updateWith(stub.article.getFirst());
		
		for (ColumnValue cv : stub.basket)
			baskets.getByName(cv.value).updateWith(cv);
	}
	
	// Possible to set by a map of parameter of a kind "property1_condition2_unit3 = K" or "property1 = LogPow
	// Comes from a HttpRequest, but can also come from somewhere else
	public void remapFromParameterMap(Map<String, String[]> parameters)
	{
		for (String key : parameters.keySet()) 
		{
			String[] value = parameters.get(key);
			RemappedValue rv = getRemappedValue(key); 
			if (rv != null)
				rv.name = value[0];
		}
	}
	
	private RemappedValue getRemappedValue(String key) //Make hierarchical, put same methods to child schemas?
	{
		String[] pieces = key.split("_");
		if (pieces[0].startsWith("property"))
		{
			PropertyRemapping pr = properties.get(Integer.valueOf(pieces[0].replace("property", "")) - 1);
			if (pieces.length == 1)
				return pr;
			if (pieces[1].startsWith("condition"))
			{
				ConditionRemapping cr = pr.conditions.get(Integer.valueOf(pieces[1].replace("condition", "")) - 1);
				if (pieces.length == 2)
					return cr;
				if (pieces[2].startsWith("unit")) 
					return cr.units.get(Integer.valueOf(pieces[2].replace("unit", "")) - 1);
				else
				if (pieces[2].startsWith("option"))
					return cr.options.get(Integer.valueOf(pieces[2].replace("option", "")) - 1);
			} else
			if (pieces[1].startsWith("unit")) 
				return pr.units.get(Integer.valueOf(pieces[1].replace("unit", "")) - 1);
			else
			if (pieces[1].startsWith("option"))
				return pr.options.get(Integer.valueOf(pieces[1].replace("option", "")) - 1); 
				
		} else
		if (pieces[0].startsWith("article"))
			return articles.get(Integer.valueOf(pieces[0].replace("article", "")) - 1);
		else
		if (pieces[0].startsWith("basket"))
			return baskets.get(Integer.valueOf(pieces[0].replace("basket", "")) - 1);
		return null;
	}
	
	public String toString()
	{
		return "properties: "+properties+"; articles: "+articles+"; baskets: "+baskets;
	}
}