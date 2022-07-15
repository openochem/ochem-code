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

import java.util.ArrayList;


public class SchemaList<T extends RemappedValue> extends ArrayList<T>
{
	private static final long serialVersionUID = 1L;
	private Class<T> genericClass;
	public T getByName(String name)
	{
		try
		{
			for (T s : this) 
				if (s.originalName.equals(name))
					return s;
			T t = genericClass.newInstance();
			t.setNames(name);
			add(t);
			return t;
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return null;
		}
	}
	public SchemaList<T> setGenericClass(Class<T> genericClass)
	{
		this.genericClass = genericClass;
		return this;
	}
}