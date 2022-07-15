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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;

/**
 * A flexible abstraction of resources of two types: classpath resources and filesystem resources.
 * The classpath resources are distinguished by the "classpath:" prefix, e.g. "classpath:models/model.xml"
 * 
 * @author midnighter
 *
 */
public class ResourceLoader {
	public static InputStream getInputStreamByPath(String path) throws FileNotFoundException
	{
		InputStream is = null;
		if (path.startsWith("classpath:"))
			is = OCHEMUtils.class.getClassLoader().getResourceAsStream(path.substring(path.startsWith("classpath:/") ? 11 : 10));
		else
			is = new FileInputStream(new File(path));
		return is;
	}

	public static Object unmarshallByPath(String path, JAXBContext jaxbContext) throws IOException, JAXBException
	{
		Object res = null;
		InputStream is = getInputStreamByPath(path);
		if (is != null)
		{
			res = jaxbContext.createUnmarshaller().unmarshal(is);
			is.close();
		}

		return res;
	}

	public static boolean resourceExists(String path) {
		if (path.startsWith("classpath:"))
			return OCHEMUtils.class.getClassLoader().getResource(path.substring(path.startsWith("classpath:/") ? 11 : 10)) != null;
		else
			return new File(path).exists();
	}

	/*
    private static String getClasspathRef(String fullPath) {
    	return fullPath.substring(fullPath.startsWith("classpath:/") ? 11 : 10);
    }
	 */
}
