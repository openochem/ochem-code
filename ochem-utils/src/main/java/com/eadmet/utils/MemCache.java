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

import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * An abstract in-memory cache.
 * Basically, an extension of a hash map to support smart clean-ups taking into account last access time
 * 
 * @author midnighter
 */
public class MemCache<T1, T2>
{
	private Map<T1, T2> map = new HashMap<T1, T2>();
	private Map<T1, Long> lastAccess = new HashMap<T1, Long>();
	
	public T2 get(T1 key) {
		T2 res = map.get(key);
		if (res != null)
			touch(key);
		return res;
	}
	
	public void put(T1 key, T2 value) {
		map.put(key, value);
		touch(key);
		
		logger.info("Adding entry " + key + " to memory cache");
	}
	
	public void cleanup(int maxElems) {
		List<Entry<T1, Long>> entries = new ArrayList<Map.Entry<T1,Long>>();
		entries.addAll(lastAccess.entrySet());
		Collections.sort(entries, new Comparator<Entry<T1, Long>>()
		{
			@Override
			public int compare(Entry<T1, Long> arg0, Entry<T1, Long> arg1)
			{
				return arg0.getValue().compareTo(arg1.getValue());
			}
		});
		
		int remove = entries.size() - maxElems;
		for (int i = 0; i < remove; i++)
		{
			T1 key = entries.get(i).getKey();
			logger.info("Removing entry " + key + " from memory cache");
			lastAccess.remove(key);
			map.remove(key);
		}
	}
	
	private void touch(T1 key) {
		lastAccess.put(key, Calendar.getInstance().getTimeInMillis());	
	}
	
	private static final Logger logger = LoggerFactory.getLogger(MemCache.class);
}
