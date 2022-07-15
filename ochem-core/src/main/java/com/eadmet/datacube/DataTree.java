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

package com.eadmet.datacube;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlAnyElement;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

/**
 * A tree of measurable data grouped by arbitrary analytical levels.
 * Can be obtained from DataCube and marshalled to XML for UI
 * 
 * @author midnighter
 *
 */
@XmlRootElement
@SuppressWarnings({"rawtypes","unchecked"})
public class DataTree<T>
{
	@XmlElement
	public List<Subkonto> subkontos;

	@XmlTransient
	public Object currentGrouping;

	@XmlElement
	public T metrics;

	@XmlElement(name = "child")
	public List<DataTree> children = new ArrayList<DataTree>();

	public void walk(Callback c) {
		c.apply(this);
		for (DataTree<T> child : children)
			child.walk(c);
	}

	@XmlAttribute
	public String getId() throws IllegalArgumentException, SecurityException, IllegalAccessException, NoSuchFieldException {
		if (currentGrouping instanceof Boolean)
			return "" + currentGrouping;
		if (currentGrouping != null)
			return "" + currentGrouping.getClass().getField("id").get(currentGrouping);
		else
			return null;
	}

	@XmlElement
	public String getTitle() {
		if (currentGrouping == null)
			return "No value";
		return "" + currentGrouping;
	}

	public static interface Callback<T> {
		public void apply(T object);
	}

	@XmlAnyElement
	private Object getCurrentGrouping() {
		if (currentGrouping instanceof Boolean)
			return null;
		return currentGrouping;
	}
}
