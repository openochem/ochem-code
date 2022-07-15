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

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.reflections.Reflections;

import com.eadmet.utils.config.ReflectionsSingleton;


public class JAXBContextFactory {
	private static final Logger logger = LogManager.getLogger(JAXBContextFactory.class);
	public static Class<?>[] getClasses()
	{
		return getClasses("");
	}

	public static Class<?>[] getClasses(String packages)
	{
		return getClasses(packages.split(":"));
	}

	public static Class<?>[] getClasses(String[] packages)
	{
		long timer = System.nanoTime();
		Reflections r = ReflectionsSingleton.get();
		Set<Class<?>> classSet = r.getTypesAnnotatedWith(XmlRootElement.class);

		Set<Field> fields = r.getFieldsAnnotatedWith(XmlElement.class);
		for (Field field : fields)
			classSet.add(field.getDeclaringClass());

		Set<Method> methods = r.getMethodsAnnotatedWith(XmlElement.class);
		for (Method method : methods)
			classSet.add(method.getDeclaringClass());

		List<Class<?>> classList = new ArrayList<Class<?>>();
		classList.addAll(classSet);

		int totalSize = classList.size();
		int i = 0;
		if (packages[0] != "")
			while (i < classList.size())
			{
				String packageName = classList.get(i).getPackage().getName();
				boolean match = false;
				for (String p : packages) 
					if (packageName.startsWith(p))
					{
						match = true;
						break;
					}

				if (match)
					i++;
				else
					classList.remove(i);
			}
		logger.info(String.format("JAXB Annotated classes: %d filtered classes out of %d total classes in %dms", classList.size(), totalSize, (System.nanoTime() - timer)/1000/1000));
		return classList.toArray(new Class<?>[1]);
	}

	public static JAXBContext get(String packages) throws JAXBException
	{
		return JAXBContext.newInstance(getClasses(packages));
	}
}
