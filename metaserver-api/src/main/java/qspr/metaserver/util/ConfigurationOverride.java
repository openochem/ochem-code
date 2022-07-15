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

package qspr.metaserver.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class ConfigurationOverride
{
	private static transient final Logger logger = LogManager.getLogger(ConfigurationOverride.class);
	
	public Map<String, Map<String, String>> properties = new HashMap<String, Map<String, String>>();
	
	public ConfigurationOverride(String fileName) throws IOException
	{
		if (!new File(fileName).exists())
		{
			logger.info("WARNING: " + fileName + " does not exist. Using the default configuration");
			return;
		}
		BufferedReader reader = new BufferedReader(new FileReader(new File(fileName)));
		
		String st;
		while ((st = reader.readLine()) != null)
		{
			st = st.trim();
			if ("".equals(st) || st.startsWith("#"))
				continue;
				
			String[] parts = st.split("=", 2);
			if (parts.length == 2)
				setProperty(parts[0].trim(), parts[1].trim());
			else
				logger.info("WARNING: Invalid configuration option " + st);
		}
		
		reader.close();
	}
	
	public ConfigurationOverride merge(ConfigurationOverride src)
	{
		for (String property : src.properties.keySet())
		{
			Map<String, String> remoteOptions = src.properties.get(property);
			Map<String, String> localOptions = properties.get(property);
			if (localOptions != null) //The property already exists, need smarter merge
			{
				for (String option : remoteOptions.keySet())
				{
					String value = remoteOptions.get(option);
					localOptions.put(option, value); //Just override what was there
				}
			} else
				properties.put(property, remoteOptions);
		}
		return this;
	}
	
	public void setProperty(String key, String value)
	{
		String parts[] = key.split("\\.", 2);
		getProperties(parts[0]).put(parts[1], value);
		logger.info("Setting configuration option: " + key + " = " + value);
	}
	
	public Map<String, String> getProperties(String parent)
	{
		Map<String, String> props = properties.get(parent);
		if (props == null)
			properties.put(parent, props = new HashMap<String, String>());
		
		return props;
	}
	
}
