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

package com.eadmet.utils.config;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlValue;
import javax.xml.bind.annotation.adapters.XmlAdapter;

interface ConfigurationItem
{
	public String getName();
}

@XmlRootElement(name="configuration")
public class ConfigurationSet implements ConfigurationItem
{
	@XmlAttribute
	public String name;
	
	@XmlAttribute
	public String comment;
	
    @XmlElement
    List<ConfigurationValue> value = new ArrayList<ConfigurationValue>();
    
    @XmlElement
    List<ConfigurationSet> set = new ArrayList<ConfigurationSet>();	
	
    
    public void addValue(ConfigurationValue value)
    {
    	this.value.add(value);
    }
    
    public void addSet(ConfigurationSet set)
    {
    	this.set.add(set);
    }

    public ConfigurationValue getValue(String name)
    {
    	for (ConfigurationValue value : this.value)
    		if (name.equals(value.name))
    			return value;
    	return null;
    }
    
    public ConfigurationSet getSet(String name)
    {
    	for (ConfigurationSet set : this.set)
    		if (name.equals(set.name))
    			return set;
    	return null;
    }
    
	public ConfigurationSet()
	{
		
	}
	
	public ConfigurationSet(String name, String comment)
	{
		this.name = name;
		this.comment = comment;
	}
	
	public ConfigurationSet(String name)
	{
		this.name = name;
	}

	public String getName() 
	{
		return name;
	}
	
	public void mergeWith(ConfigurationSet target)
	{
		for (ConfigurationValue value : target.value) 
		{
			ConfigurationValue local = getValue(value.name);
			if (local == null)
				addValue(value);
			else
				local.value = value.value;
		}
		
		for (ConfigurationSet set : target.set)
		{
			ConfigurationSet local = getSet(set.name);
			if (local == null)
				addSet(set);
			else
				local.mergeWith(set);
		}
	}
	
	public Map<String, String> getValuesMap()
	{
		Map<String, String> map = new HashMap<String, String>();
		for (ConfigurationValue value : this.value)
			map.put(value.name, value.value);
		return map;
	}
}

@XmlRootElement
class ConfigurationValue implements ConfigurationItem
{
	
	@XmlAttribute
	public String name;

	@XmlAttribute
	public String comment;
	
	@XmlValue
	public String value;
	
	public ConfigurationValue()
	{
		
	}
	
	public ConfigurationValue(String name, String value, String comment)
	{
		this.name = name;
		this.value = value;
		this.comment = comment;
	}
	
	public ConfigurationValue(String name, String value)
	{
		this.name = name;
		this.value = value;
	}
	
	public String getName()
	{
		return name;
	}
}

class ConfigurationValueMapAdapter extends XmlAdapter<ConfigurationValue[], Map<String, ConfigurationValue>>
{
    public ConfigurationValue[] marshal(Map<String, ConfigurationValue> arg) throws Exception 
    {
    	ConfigurationValue[] mapElements = new ConfigurationValue[arg.size()];
        int i=0;
    	for (ConfigurationValue item : arg.values())
    	{
    		mapElements[i] = item;
    		i++;
    	}

        return mapElements;
    }

    public Map<String, ConfigurationValue> unmarshal(ConfigurationValue[] arg) throws Exception
    {
        Map<String, ConfigurationValue> r = new HashMap<String, ConfigurationValue>();
        
        for (ConfigurationValue item : arg)
            r.put(item.getName(), item);
        		
        return r;
    }
}

class ConfigurationSetMapAdapter extends XmlAdapter<ConfigurationSet[], Map<String, ConfigurationSet>>
{
    public ConfigurationSet[] marshal(Map<String, ConfigurationSet> arg) throws Exception 
    {
    	ConfigurationSet[] mapElements = new ConfigurationSet[arg.size()];
        int i=0;
    	for (ConfigurationSet item : arg.values())
    	{
    		mapElements[i] = item;
    		i++;
    	}

        return mapElements;
    }

    public Map<String, ConfigurationSet> unmarshal(ConfigurationSet[] arg) throws Exception
    {
        Map<String, ConfigurationSet> r = new HashMap<String, ConfigurationSet>();
        
        for (ConfigurationSet item : arg)
            r.put(item.getName(), item);
        		
        return r;
    }
}
