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

import java.util.Iterator;
import java.util.Properties;
import java.util.Set;

import com.eadmet.exceptions.UserFriendlyException;

public final class OSType
{
	private static String operatingSystem = null;
	private static String pathSeparator = null;
	private static String fileSeparator = null;
	private static String endLine = null;

	private static String PREFIXMAC = "/opt/local/bin/python";
	private static String PREFIXLIN = "/opt/conda/bin/python";


	public static String getPython(int version) {
		switch(version) {
		case 39: return (OSType.isMac()?PREFIXMAC:PREFIXLIN) + version/10.;
		}
		throw new UserFriendlyException("Required python version is not found version: " + version);
	}

	public static String getOsName()
	{
		if (operatingSystem == null)
			operatingSystem = System.getProperty("os.name");
		return operatingSystem;
	}

	// Why do we need this? If possible, everything OS-specific should be encapsulated here
	public static boolean isWindows()
	{
		return getOsName().startsWith("Windows");
	}

	// Why do we need this? If possible, everything OS-specific should be encapsulated here
	public static boolean isMac()
	{
		return getOsName().startsWith("Mac");
	}

	public static String endLine()
	{
		if (endLine != null)
			return endLine;
		return endLine = System.getProperty("line.separator");
	}
	public static String pathSeparator()
	{
		if (pathSeparator != null)
			return pathSeparator;
		return pathSeparator = System.getProperty("path.separator");
	}

	public static String fileSeparator()
	{
		if (fileSeparator != null)
			return fileSeparator;
		return fileSeparator = System.getProperty("file.separator");
	}

	/**
	 * Glue a path using the OS-specific separator
	 * 
	 * @author midnighter
	 * @param pieces
	 * @return
	 */
	public static String getPath(String... pieces)
	{
		String result = "";
		for (int i = 0; i < pieces.length; i++)
		{
			result += pieces[i];
			if (i < pieces.length - 1)
				result += fileSeparator();
		}

		return result;
	}

	public static String convertPath(String path) {
		if (OSType.isWindows())
			return path.replaceAll("/", "\\\\");
		return path;
	}

	public static String getHome(){
		return (getOsName().startsWith("Mac")?"/Users/":"/home/")+System.getProperty("user.name")+"/";
	}

	public static void main(String[] args)
	{
		Properties prop = System.getProperties();
		Set<String> a = prop.stringPropertyNames();
		Iterator<String> keys = a.iterator();
		while (keys.hasNext())
		{
			String key = keys.next();
			String value = System.getProperty(key);
			System.out.println(key + "=" + value);
		}
	}

	public static boolean isAarch64() {
		return System.getProperty("os.arch").contains("aarch64");
	}

}
