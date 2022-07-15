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

package qspr.business;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import com.eadmet.utils.OCHEMUtils;


@XmlRootElement(name = "filters")
public class WebFilters 
{
	@XmlTransient
	private Map<String, Filter> internal = new HashMap<String, Filter>();
	
	@XmlElement
	public List<Filter> filter;
	
	public void addFilter(String name, String value, String title)
	{		
		if (!has(name))
		{
			addFilterOverride(name, value, title);
		}
	}

	public void addFilterOverride(String name, String value, String title)
	{
		if (filter == null)
			filter = new ArrayList<Filter>();
		
		Filter curfilter;
		
		if (has(name))
		{
			curfilter = internal.get(name);
			curfilter.value = value;
		} else
		{
			curfilter = new Filter();
			curfilter.name = name;
			curfilter.title = title;
			curfilter.value = value;
			filter.add(curfilter);
			internal.put(name, curfilter);
		}
	}
	
	public boolean has(String name)
	{
		return internal.containsKey(name) && !"".equals(internal.get(name).value);
	}
	
	public void unset(String name)
	{
		if (has(name))
		{
			//filter.remove(get(name));
			filter.remove(internal.get(name)); //TODO Check that my correction is correct one
			internal.remove(name);
		}
	}
	
	public boolean notEmpty(String name)
	{
		return has(name) && !get(name).equals(""); 
	}
	
	public String get(String name)
	{
		if (has(name))
			return internal.get(name).value;
		else 
			return null;
	}
	
	public Integer getInteger(String name)
	{
		String names = get(name);
		names = names.trim();
		return Integer.valueOf(names);
	}
	
	public Long getLong(String name)
	{
		String names = get(name);
		names = names.trim();
		return Long.valueOf(names);
	}
	
	public Long[] getLongArray(String name)
	{
		String[] parts = get(name).split(",");
		Long[] res = new Long[parts.length];
		for (int i = 0; i < res.length; i++)
			res[i] = Long.valueOf(parts[i].trim());
		return res;
	}
	
	public Set<String> getFilterSet()
	{
		return internal.keySet();
	}
	
	public boolean in(String name, String[] values)
	{
		for (String value : values) {
			if (value.equals(get(name)))
				return true;
		}
		
		return false;
	}
	
	public String getHash()
	{
		Collections.sort(filter);
		StringBuffer buff = new StringBuffer();
		for (Filter f : filter)
			if (!"pagenum".equals(f.name))
			{
				buff.append(f.name);
				buff.append(f.value);
			}
		
		return OCHEMUtils.getMD5(buff.toString());
	}
}

class Filter implements Comparable<Filter>
{
	@XmlAttribute
	String name;
	
	@XmlAttribute
	
	String value;
	
	@XmlAttribute
	String title;
	
	@XmlElement(name = "option")
	List<Filter> options;

	public int compareTo(Filter compareTo)
	{
		return name.compareTo(compareTo.name);
	}
}
