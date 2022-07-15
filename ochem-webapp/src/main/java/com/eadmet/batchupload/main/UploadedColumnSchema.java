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
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement
public class UploadedColumnSchema 
{
	@XmlAttribute
	public Boolean ignore = false;

	@XmlAttribute
	public Boolean hidden = false;

	@XmlAttribute
	public String name;

	@XmlAttribute
	public String originalName;

	@XmlAttribute
	public String subscriptName;


	@XmlAttribute
	public Integer realNumber;

	@XmlAttribute
	public Integer remappedNumber;

	@XmlAttribute
	public String defaultValue;

	@XmlElement
	public  List<String> sampleValues = new ArrayList<String>();

	@XmlAttribute
	public UploadedColumnType type;

	public UploadedColumnSchema()
	{
	}

	public UploadedColumnSchema(String columnName, Integer columnNumber)
	{
		this.name = this.originalName = columnName;
		this.realNumber = columnNumber;
		this.type = UploadedColumnType.undefined;
	}

	public UploadedColumnSchema setHidden(boolean hidden)
	{
		this.hidden = hidden;
		return this;
	}

	public String toString()
	{
		return realNumber+":"+remappedNumber+":"+name+":"+originalName+":"+type+":"+ignore;
	}


	public enum UploadedColumnType
	{
		//Try.... remove predicate as a separate column
		undefined,
		molecule("molecule", "smile", "smiles", "sdf", "moleculeid", "mixture"),
		predicate("predicate"),
		article("article", "articleid", "pubmedid"),
		property("property"),condition("condition"),
		property_and_value, condition_and_value, //Legacy column types with properties as headers
		value("value","min_value"), accuracy("accuracy"),  max_value("max_value","interval"), 
		unit("unit"), 
		page("page"), line("line"), table("table"), n("n"), hidden("hidden"),
		name("name"), casrn("casrn","cas"), externalid("externalid","external_id"), 
		recordid("recordid"), 
		basket("basket","moleculeset"), comment("comment","comments"), evidence("evidence");

		String[] names;
		static Map<String, UploadedColumnType> nameMap;

		static UploadedColumnType fromString(String v)
		{
			prepareMap();
			v = v.toLowerCase();
			if (nameMap.get(v) != null)
				return nameMap.get(v);
			else
				return undefined;
		}

		private static void prepareMap()
		{
			if (nameMap != null)
				return;

			nameMap = new HashMap<String, UploadedColumnType>();
			for (UploadedColumnType type : UploadedColumnType.values())
				for (String name : type.names) 
					nameMap.put(name, type);
		}

		private UploadedColumnType(String... names) //
		{
			this.names = names;
		}

		public static List<String> getKnownColumns()
		{
			List<String> l = new ArrayList<String>();
			for (UploadedColumnType type : UploadedColumnType.values())
				if (type.names.length > 0)
					l.add(type.names[0]);
			l.add("mixture"); // additional names
			Collections.sort(l);
			return l;
		}
	}
}

