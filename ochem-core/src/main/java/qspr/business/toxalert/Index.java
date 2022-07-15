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

package qspr.business.toxalert;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * A very simple index. Its a helper class for fast display of the filtered structural alerts
 * 
 * @author midnighter
 *
 * @param <T1>
 * @param <T2>
 */

// Midnighter on Mar 12, 2012
public class Index<T1, T2> {
	public List<T1> index = new ArrayList<T1>();
	public List<List<T2>> values = new ArrayList<List<T2>>();
	
	// This maps are for speed. Can be refactored. 
	// The question here is how to combine indexing features of a Map/Set with a positional access of List. Needs to be investigated. So far, use both Set/Map and List
	public List<Set<T2>> valuesSet = new ArrayList<Set<T2>>(); 
	public Map<T1, Integer> indexPositionHash = new HashMap<T1, Integer>();

	public int addValue(T1 indexValue, T2 value) {
		Integer i = indexPositionHash.get(indexValue); 
		if (i == null) {
			index.add(indexValue);
			values.add(new ArrayList<T2>());
			valuesSet.add(new HashSet<T2>());
			i = index.size() - 1;
			indexPositionHash.put(indexValue, i);
		}
		if (!valuesSet.get(i).contains(value))
		{
			values.get(i).add(value);
			valuesSet.get(i).add(value);
		}
		return i;
	}
	
	/**
	 * Get an inverse index
	 * @return
	 */
	public Index<T2, T1> getInverse()
	{
		Index<T2, T1> inverse = new Index<T2, T1>();
		
		for (int i = 0; i < index.size(); i++)
			for (int k = 0; k < values.get(i).size(); k++)
				inverse.addValue(values.get(i).get(k), index.get(i));
		
		return inverse;
	}
}