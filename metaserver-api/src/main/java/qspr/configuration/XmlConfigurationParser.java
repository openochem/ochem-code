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

package qspr.configuration;

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.Reader;
import java.io.Writer;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

import com.eadmet.utils.config.ConfigurationSet;

public class XmlConfigurationParser 
{
	@SuppressWarnings("rawtypes")
	static Class[] classes = {ConfigurationSet.class};
	static JAXBContext jaxbContext;
	static 
	{
		try
		{
			jaxbContext = JAXBContext.newInstance(classes);
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	public static ConfigurationSet parseConfiguration(Reader r) throws Exception
	{
		Unmarshaller um = jaxbContext.createUnmarshaller();
		ConfigurationSet set = (ConfigurationSet)um.unmarshal(r);
		return set;
	}

	public static ConfigurationSet parseConfiguration(String[] fileNames) throws Exception
	{
		ConfigurationSet set = new ConfigurationSet();
		for (String fileName : fileNames) 
		{
			File f = new File(fileName); 
			if (!f.exists())
				continue;
			ConfigurationSet cs = (ConfigurationSet)parseConfiguration(new FileReader(f));
			set.mergeWith(cs);
		}
		return set;
	}

	public static void exportConfiguration(ConfigurationSet set, Writer w) throws Exception
	{
		Marshaller m = jaxbContext.createMarshaller();
		m.setProperty("jaxb.formatted.output",Boolean.TRUE);
		m.marshal(set, w);
	}

	public static void exportConfiguration(ConfigurationSet set, String fileName) throws Exception
	{
		FileWriter fw = new FileWriter(fileName);
		exportConfiguration(set, fw);
		fw.close();
	}
}

