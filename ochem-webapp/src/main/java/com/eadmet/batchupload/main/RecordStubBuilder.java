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

package com.eadmet.batchupload.main;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.eadmet.batchupload.main.RecordStub.ColumnValue;
import com.eadmet.batchupload.main.RecordStub.NameValueStub;
import com.eadmet.batchupload.main.RecordStub.ValueStub;
import com.eadmet.exceptions.UserFriendlyException;

import qspr.dao.Repository;
import qspr.entities.Property;

/*
 * A class that converts a UpoadedRow obtained from the parser to a RecordStub (or several, if the row contains multiple records)
 * The processColumnValue is a huge hard-coded domain-specific and grammar-specific function that fills the appropriate stub fields
 * 
 * The hard-coded rules are: each "property" column denotes a new record, molecule and article colums are common for all records
 * The "value", "unit", etc. columns belong to the nearest (appeared) property or condition
 */
public class RecordStubBuilder 
{

	Map<String,String> defaultUnits = new HashMap<String, String>();

	List<RecordStub> stubs = null; // List of "records" produced from the row
	RecordStub r = null; // Current record to fill with column values

	List<ValueStub> uvs = new ArrayList<ValueStub>(); //Unassigned value stubs

	ValueStub vs; //Current value stub to put child fields
	NameValueStub nvs; //Current name value stub (condition) to put child nodes

	public List<RecordStub> getRowStubs(UploadedRow row)
	{
		stubs = new ArrayList<RecordStub>();
		r = new RecordStub();
		uvs = new ArrayList<ValueStub>();
		vs = new ValueStub();
		uvs.add(vs);

		stubs.add(r);
		for (UploadedColumnSchema column : row.parentSheet.getColumns()) 
		{
			if (column.ignore)
				continue;

			String value = null;
			if (column.realNumber != null)
				value = row.values == null || row.values.size() <= column.realNumber? null : row.values.get(column.realNumber);
			else
				value = column.defaultValue; // this can be only Dummy property; thus default value will be used for it

			if (value == null){
				value = column.defaultValue;
				continue;
			}

			processColumnValue(column, value);
		}
		setDefaults(row);
		return stubs;
	}

	private String getDefaultUnit(String name) {
		if(defaultUnits.containsKey(name))return defaultUnits.get(name);
		Property prop = Repository.property.getProperty(name, false);
		defaultUnits.put(name,prop == null ? "default": prop.defaultUnit.getName());
		return defaultUnits.get(name);		
	}

	private void setDefaults(UploadedRow row)
	{
		for (RecordStub stub : stubs) 
		{
			if (stub.property.value == null || stub.property.value.size() == 0)
				continue;

			if (stub.n.size() == 0)
				stub.n.add(new ColumnValue(null, "AUTO_"+(row.number+1)));

			if (stub.property.value.getFirst().unit.size() == 0)
				stub.property.value.getFirst().unit.add(new ColumnValue(null, getDefaultUnit(stub.property.name.value)));

			for (NameValueStub nvs : stub.property.conditions)
			{
				if(nvs.value == null || nvs.value.size() == 0)
					throw new UserFriendlyException("\n\n\nDear User. Check your data format.\n"
							+ "In the technical format EACH specified condition should have A VALUE.\n"
							+ "If value is not available, do not provide the respective condition name.\n"
							+ "This is not a bug and it will not be fixed in the nearest future.\n\n\n");

				if (nvs.value.getFirst().unit.size() == 0)
					nvs.value.getFirst().unit.add(new ColumnValue(null, getDefaultUnit(stub.property.name.value)));	
			}

			if (stub.article.size() == 0)
				stub.article.add(new ColumnValue(null, "unpublished"));

			boolean addDefault = true;
			for (ColumnValue cv : stub.basket) 
				if (cv.value.equals("default"))
					addDefault = false;
			if (addDefault)
				stub.basket.add(new ColumnValue(null, "default"));
		}
	}

	private void processColumnValue(UploadedColumnSchema column, String value)
	{

		switch (column.type)
		{
		case accuracy:
			vs.accuracy.add(new ColumnValue(column, value));
			break;
		case predicate:
			vs.predicate.add(new ColumnValue(column, value));
			break;
		case article:
			r.article.add(new ColumnValue(column, value));
			break;
		case basket:
			r.basket.add(new ColumnValue(column, value));
			break;
		case comment:
			r.comments.add(new ColumnValue(column, value));
			break;
		case condition:
			if (vs.value != null)
			{
				vs = new ValueStub();
				uvs.add(vs);
			}				
			if (r.property.conditions.size() == 0 || r.property.conditions.getLast().name != null)
				r.property.conditions.add(new NameValueStub());
			r.property.conditions.getLast().name = new ColumnValue(column, value);
			r.property.conditions.getLast().value.addAll(uvs);
			uvs.clear();
			nvs = r.property.conditions.getLast();
			break;
		case condition_and_value:
			if (vs.value != null)
			{
				vs = new ValueStub();
				uvs.add(vs);
			}
			vs.value = new ColumnValue(column, value);

			if (r.property.conditions.size() == 0 || r.property.conditions.getLast().name != null)
				r.property.conditions.add(new NameValueStub());
			r.property.conditions.getLast().name = new ColumnValue(column, column.name);

			r.property.conditions.getLast().value.addAll(uvs);
			uvs.clear();
			nvs = r.property.conditions.getLast();
			break;
		case evidence:
			r.evidence.add(new ColumnValue(column, value));
			break;
		case externalid:
			r.externalid.add(new ColumnValue(column, value));
			break;
		case hidden:
			r.hidden.add(new ColumnValue(column, value));
			break;
		case line:
			r.line.add(new ColumnValue(column, value));
			break;
		case max_value:
			vs.interval.add(new ColumnValue(column, value));
			break;
		case molecule:
			r.molecules.add(new ColumnValue(column, value));
			break;
		case n:
			r.n.add(new ColumnValue(column, value));
			break;
		case name:
		case casrn: // casrn the same as name!
			String names[] = value.split(";");
			for(String n:names)
				r.name.add(new ColumnValue(column, n.trim()));
			break;
		case page:
			r.page.add(new ColumnValue(column, value));
			break;
		case property:
			if (vs.value != null)
			{
				vs = new ValueStub();
				uvs.add(vs);
			}
			if (r.property.name != null)
			{
				r = new RecordStub(r);
				stubs.add(r);
			}
			r.property.name = new ColumnValue(column, value);
			r.property.value.addAll(uvs);
			uvs.clear();
			nvs = r.property;
			break;
		case property_and_value:

			//FIXME Fix to work with to values; really, something better is required!
			if(value.indexOf(" to ")!=-1){
				value = value.replaceAll(" to ", " - ");
			}

			if (vs.value != null)
			{
				vs = new ValueStub();
				uvs.add(vs);
			}
			vs.value = new ColumnValue(column, value);
			if (r.property.name != null)
			{
				r = new RecordStub(r);
				stubs.add(r);
			}
			r.property.name = new ColumnValue(column, column.name);

			r.property.value.addAll(uvs);
			uvs.clear();
			nvs = r.property;
			break;
		case recordid:
			break;
		case table:
			r.table.add(new ColumnValue(column, value));
			break;
		case undefined:
			break;
		case unit:
			vs.unit.add(new ColumnValue(column, value));
			break;
		case value:
			if (vs.value != null)
			{
				vs = new ValueStub();
				if (nvs == null)
					uvs.add(vs);
				else
					nvs.value.add(vs);
			}
			vs.value = new ColumnValue(column, value);
			break;
		default:
			break;
		}
	}
}
