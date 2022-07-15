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

package com.eadmet.mmpa;

import java.io.Serializable;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import com.eadmet.utils.PageInfo;

public class TransformationFilter implements Serializable
{
	public Long minCount;
	
	public Long exactCount;
	
	public Long fragId;
	
	public boolean splitSets;
	/**
	 * Filter the transformations matching particular molecule.
	 * This is an invariant identifier (M-ID, "Mapping2" ID)
	 */
	public Integer moleculeId;
	
	/**
	 * Explicit list of transformation IDs
	 */
	public List<Long> ids;
	
	/**
	 * A map of desired effects for a given property in a given transformation set
	 * An array of booleans stands for {decreasing, no effect, increasing} effect flags
	 */
	public Map<Long, Map<Long, boolean[]>> propertyEffect = new HashMap<Long, Map<Long, boolean[]>>();
	
	PageInfo pageInfo;
	
	public boolean[] getEffects(Long setId, Long propertyId) 
	{
		Map<Long, boolean[]> m = propertyEffect.get(setId);
		if (m == null)
			propertyEffect.put(setId, m = new HashMap<Long, boolean[]>());
		
		boolean[] effects = m.get(propertyId);
		
		if (effects == null)
			m.put(propertyId, effects = new boolean[3]);
		
		return effects;
	}
	
	/**
	 * Get the desired effects for the given property in the first suitable set
	 */
	public boolean[] getEffects(Long propertyId) {
		for (Long setId : propertyEffect.keySet())
			if (propertyEffect.get(setId).containsKey(propertyId))
				return propertyEffect.get(setId).get(propertyId);
		return null;
	}
	
	public void clean() 
	{
		Iterator<Entry<Long, Map<Long, boolean[]>>> i = propertyEffect.entrySet().iterator();
		while (i.hasNext())
		{
			Entry<Long, Map<Long, boolean[]>> e = i.next();
			if (e.getValue().isEmpty())
			{
				i.remove();
				continue;
			}
			Iterator<Entry<Long, boolean[]>> iter = e.getValue().entrySet().iterator();
			while (iter.hasNext())
			{
				Entry<Long, boolean[]> entry = iter.next();
				if (entry.getValue()[0] && entry.getValue()[1] && entry.getValue()[2])
					iter.remove();
			}
		}
	}
	
	public Set<Long> getPropertyIds() {
		Set<Long> ids = new HashSet<Long>();
		for (Long set : propertyEffect.keySet())
			for (Long id : propertyEffect.get(set).keySet())
				ids.add(id);
		
		return ids;
	}
	
	public List<Long> getSetsByProperty(Long propertyId) {
		List<Long> sets = new ArrayList<Long>();
		for (Long setId : propertyEffect.keySet())
			if (propertyEffect.get(setId).containsKey(propertyId))
				sets.add(setId);
		return sets;
	}
	
	public String toString() {
		StringWriter sw = new StringWriter();
		sw.append("" + getPropertyIds().size() + " propertie(s)");
		if (moleculeId != null)
			sw.append(", molecule M" + moleculeId);
		
		return sw.toString();
	}
	
	private static final long serialVersionUID = 1L;
}