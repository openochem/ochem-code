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

package com.eadmet.utils;

import java.util.HashMap;
import java.util.Map;

/**
 * A helper utility for "maps of maps"
 * 
 * @author midnighter
 *
 * @param <T1> first key type
 * @param <T2> second key type
 * @param <T3> value type
 */
public class DoubleMap<T1, T2, T3> extends HashMap<T1, Map<T2, T3>>
{
	/**
	 * Add a value for the two provided keys
	 */
	public void put(T1 key1, T2 key2, T3 val) {
		Map<T2, T3> vals = get(key1);
		if (vals == null)
			put(key1, vals = new HashMap<T2, T3>());
		vals.put(key2, val);
	}
	
	public T3 get(T1 key1, T2 key2) {
		Map<T2, T3> vals = get(key1);
		if (vals == null)
			return null;
		return vals.get(key2);
	}
	
	public DoubleMap<T1, T2, T3> clone() {
		DoubleMap<T1, T2, T3> clone = new DoubleMap<T1, T2, T3>();
		for (T1 t1 : keySet())
			clone.put(t1, new HashMap<T2, T3>(get(t1)));
		
		return clone;
	}
	
	private static final long serialVersionUID = 1L;
}
