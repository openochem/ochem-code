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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.io.StringWriter;
import java.io.Writer;
import java.net.URL;
import java.util.ArrayList;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class PlainTextConfigurationParser 
{
	private static final Logger logger = LoggerFactory.getLogger(PlainTextConfigurationParser.class);
	
	static Class<?>[] classes = {ConfigurationSet.class};
	
	public ConfigurationSet confSet = new ConfigurationSet();
	
	public static ConfigurationSet parseConfiguration(BufferedReader ir, String resourceName) throws Exception
	{
		// Resolve includes 
		StringWriter sw = new StringWriter();
		BufferedWriter bsw = new BufferedWriter(sw);
		readStreamWithResolvedIncludes(ir, bsw, resourceName, new ArrayList<String>());
		bsw.close();
		
		BufferedReader r = new BufferedReader(new StringReader(sw.toString()));
		
		ConfigurationSet set = new ConfigurationSet();
		String s;
		int i = 0;
		while ((s = r.readLine()) != null)
		{
			i++;
			s = s.trim();
			if (s.startsWith("#") || s.equals(""))
				continue;
			String[] keyvalue = s.split("\\s*=\\s*", 2);
			
			if (keyvalue.length != 2)
				throw new Exception("Error in configuration on line "+i+": no key-value pair found");
			
			String[] keys = keyvalue[0].split("\\.(?=([^\\}]*\\{[^\\}]*\\})*[^\\{\\}]*$)");
			ConfigurationSet subset = set;
			for (int j=0; j<keys.length; j++)
			{
				String keyPiece = keys[j].replaceAll("[\\{\\}\\s]", "");
				
				if (j == keys.length-1)
				{	
					if (subset.getValue(keyPiece) != null)
						throw new Exception("Error in configuration on line "+i+": duplicate key "+keyvalue[0]);
					subset.addValue(new ConfigurationValue(keyPiece, keyvalue[1]));
				} else
				{
					ConfigurationSet tmp = subset.getSet(keyPiece);
					if (tmp == null)
					{
						tmp = new ConfigurationSet(keyPiece);
						subset.addSet(tmp);
					}
					subset = tmp;
				}
			}
		}
		return set;		
	}
	
	public static ConfigurationSet parseConfiguration(BufferedReader r) throws Exception
	{
		return parseConfiguration(r, null);
	}
	
	public static ConfigurationSet parseConfiguration(String[] fileNames) throws Exception
	{
		ConfigurationSet set = new ConfigurationSet();
		for (String fileName : fileNames) 
		{
			ConfigurationSet cs = parseConfiguration(fileName);
			if (cs != null)
				set.mergeWith(cs);
		}
		return set;
	}
	
	public static ConfigurationSet parseConfiguration(String fileName) throws Exception
	{
		if (fileName == null)
			return null;
		File f = new File(fileName);
		if (!f.exists())
			return null;
		
		logger.info("Loading configuration resource " + f.getAbsolutePath());
		ConfigurationSet cs = (ConfigurationSet)parseConfiguration(new BufferedReader(new FileReader(f)),f.getAbsolutePath());
		return cs;
	}
	private static void readStreamWithResolvedIncludes(BufferedReader r, BufferedWriter w, String resourceName, List<String> processedIncludes) throws IOException
	{
		if (processedIncludes.contains(resourceName))
			throw new IOException("Multiple include attempt or circular reference for file "+resourceName);
		
		processedIncludes.add(resourceName);
		String s = null;
		while ((s = r.readLine()) != null)
		{
			if (!s.startsWith("#include"))
			{
				w.write(s+"\n");
				continue;
			}
			String path = s.split("\\s",2)[1];
			if (resourceName != null && new File(resourceName).exists()) //Real file
			{
				File f = new File(resourceName);
				if (!path.startsWith("/"))
					path = f.getParentFile().getAbsolutePath()+"/"+path;
				File nf = new File(path);
				if (nf.exists())
					readStreamWithResolvedIncludes(new BufferedReader(new FileReader(path)), w, path, processedIncludes);
				else
					throw new FileNotFoundException("Included file not found "+nf.getAbsolutePath());
			} else //Assume it's a classloader resource
			{
				URL url = PlainTextConfigurationParser.class.getClassLoader().getResource(path);
				if (url != null)
					readStreamWithResolvedIncludes(new BufferedReader(new InputStreamReader(url.openStream())), w, path, processedIncludes);
				else
					throw new FileNotFoundException("Included resource not found "+resourceName);
			}
		}
	}
	

	
	private static void writeBlockComment(String comment, Writer w) throws IOException
	{
		{
			w.write("\n");
			for (int i=0; i<comment.length()+4; i++)
				w.write("#");
			w.write("\n# "+comment+" #\n");
			for (int i=0; i<comment.length()+4; i++)
				w.write("#");
			w.write("\n");
		}
	}
	
	private static void exportConfiguration(ConfigurationSet set, BufferedWriter w, String prefix) throws Exception
	{
		for (ConfigurationValue value : set.value) 
		{
			String strValue = value.value == null ? "" : value.value;
			if (value.comment != null && !value.comment.equals(""))
				w.write("\n# "+value.comment+"\n");
			
			if (value.name.contains("."))
				w.write(prefix+"{"+value.name+"} = "+strValue+"\n");
			else
				w.write(prefix+value.name+" = "+strValue+"\n");
			if (value.comment != null && !value.comment.equals(""))
				w.write("\n");
		}
		for (ConfigurationSet subset : set.set) 
		{
			if (subset.comment != null && !subset.comment.equals(""))
				writeBlockComment(subset.comment, w);
			
			String newPrefix = (subset.name.contains(".")) ? prefix+"{"+subset.name+"}." : prefix+subset.name+".";
			exportConfiguration(subset, w, newPrefix);
		}
	}
	
	public static void exportConfiguration(ConfigurationSet set, BufferedWriter w) throws Exception
	{
		exportConfiguration(set, w, "");
	}
	
	public static void exportConfiguration(ConfigurationSet set, String fileName) throws Exception
	{
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		exportConfiguration(set, bw);
		bw.close();
	}
}

