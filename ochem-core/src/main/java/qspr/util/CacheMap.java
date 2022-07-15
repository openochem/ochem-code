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

package qspr.util;

import java.util.HashMap;

public class CacheMap<T> extends HashMap<String, Object>
{
	private static final long serialVersionUID = 1L;

	public void put(String key, T value)
	{
		super.put(key, value);
	}

	public void put(String key, Exception value)
	{
		super.put(key, value);
	}

	@SuppressWarnings("unchecked")
	public T get(String key) throws Exception
	{
		Object o = super.get(key);

		if (o == null)
			return null;

		if (o instanceof Exception)
			throw (Exception)o;

		return (T)o;

	}
}