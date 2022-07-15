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
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.reflections.Reflections;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class GlobalConfigurator
{
	private static transient final Logger logger = LoggerFactory.getLogger(GlobalConfigurator.class);

	@SuppressWarnings("rawtypes")
	private static Map<Class, Class> classMap = new HashMap<Class, Class>()
	{
		private static final long serialVersionUID = 1L;
		{
			put(boolean.class, Boolean.class);
			put(byte.class, Byte.class);
			put(char.class, Character.class);
			put(double.class, Double.class);
			put(float.class, Float.class);
			put(int.class, Integer.class);
			put(long.class, Long.class);
			put(short.class, Short.class);
			put(Boolean.class, Boolean.class);
			put(Byte.class, Byte.class);
			put(Character.class, Character.class);
			put(Double.class, Double.class);
			put(Float.class, Float.class);
			put(Integer.class, Integer.class);
			put(Long.class, Long.class);
			put(Short.class, Short.class);
		}
	};

	private ConfigurationSet confSet = new ConfigurationSet();

	public void addResource(String fileName) throws Exception
	{
		addConfiguration(PlainTextConfigurationParser.parseConfiguration(fileName));
	}

	public void addResource(URL url) throws Exception
	{
		if (url != null)
		{
			InputStream is = url.openStream();
			addConfiguration(PlainTextConfigurationParser.parseConfiguration(new BufferedReader(new InputStreamReader(is))));
			is.close();
		}
	}

	public void configure() throws Exception
	{
		configure(confSet);
	}

	public void addConfiguration(ConfigurationSet cSet)
	{
		if (cSet != null)
			confSet.mergeWith(cSet);
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static void configure(ConfigurationSet projectConfig) throws Exception
	{
		Reflections reflections = ReflectionsSingleton.get(); //All our classes
		Set<Class<?>> annotated = reflections.getTypesAnnotatedWith(ConfigurableClass.class);
		for (Class<?> configurableClass : annotated)
		{
			String className = configurableClass.getAnnotation(ConfigurableClass.class).name();
			ConfigurationSet classConfig = projectConfig.getSet(className);

			if (classConfig == null)
				continue;

			for (Field field : configurableClass.getDeclaredFields())
			{
				if (!field.isAnnotationPresent(ConfigurableProperty.class))
					continue;

				field.setAccessible(true);
				String fieldName = field.getAnnotation(ConfigurableProperty.class).name();

				ConfigurationSet fieldConfigSet = classConfig.getSet(fieldName);
				ConfigurationValue fieldConfigValue = classConfig.getValue(fieldName);

				Class c = field.getType();

				try
				{
					if (c.equals(Map.class))
					{
						if (fieldConfigSet != null)
							field.set(null, fieldConfigSet.getValuesMap());
					}
					else if (c.equals(String.class))
					{
						if (fieldConfigValue != null)
							field.set(null, fieldConfigValue.value);
					}
					else if (classMap.containsKey(c))
					{
						Class converterClass = classMap.get(c);
						Method m = converterClass.getMethod("valueOf", String.class);
						if (fieldConfigValue != null)
							field.set(null, m.invoke(null, fieldConfigValue.value));
					}
					else
						throw new Exception("Unsupported value class " + c.getCanonicalName());
				} catch (Exception e)
				{
					logger.info("[WARNING] Can't configure field " + fieldName + " of class " + className);
					e.printStackTrace();
				}
			}
		}
	}

	public static void configure(String fileName) throws Exception
	{
		configure(new String[] { fileName });
	}

	public static void configure(String[] fileNames) throws Exception
	{
		ConfigurationSet projectConfig = PlainTextConfigurationParser.parseConfiguration(fileNames);
		configure(projectConfig);
	}

	public static void export(String fileName) throws Exception
	{
		BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
		export(bw);
		bw.close();
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static void export(BufferedWriter w) throws Exception
	{
		ConfigurationSet projectConfig = new ConfigurationSet();

		Reflections reflections = ReflectionsSingleton.get();
		Set<Class<?>> annotated = reflections.getTypesAnnotatedWith(ConfigurableClass.class);
		for (Class<?> configurableClass : annotated)
		{
			String className = configurableClass.getAnnotation(ConfigurableClass.class).name();
			String classComment = configurableClass.getAnnotation(ConfigurableClass.class).comment();
			ConfigurationSet classConfig = new ConfigurationSet(className, classComment);

			for (Field field : configurableClass.getDeclaredFields())
			{
				if (!field.isAnnotationPresent(ConfigurableProperty.class))
					continue;

				field.setAccessible(true);
				String fieldName = field.getAnnotation(ConfigurableProperty.class).name();
				String fieldComment = field.getAnnotation(ConfigurableProperty.class).comment();

				Class c = field.getType();

				try
				{
					if (c.equals(Map.class))
					{
						Map<String, String> value = (Map<String, String>) field.get(null);

						ConfigurationSet set = new ConfigurationSet(fieldName, fieldComment);

						if (value != null)
							for (String key : value.keySet())
								set.addValue(new ConfigurationValue(key, value.get(key)));

						classConfig.addSet(set);
					}
					else if (c.equals(String.class))
					{
						classConfig.addValue(new ConfigurationValue(fieldName, (String) field.get(null), fieldComment));
					}
					else if (classMap.containsKey(c))
					{
						Object value = field.get(null);
						if (value != null)
						{
							Method m = value.getClass().getMethod("toString", (Class[]) null);
							classConfig.addValue(new ConfigurationValue(fieldName, (String) m.invoke(value, (Object[]) null), fieldComment));
						}
						else
							classConfig.addValue(new ConfigurationValue(fieldName, null, fieldComment));
					}
					else
						throw new Exception("Unsupported type " + c.getCanonicalName());
				} catch (Exception e)
				{
					logger.info("[WARNING] Can't export field " + fieldName + " of class " + className);
					e.printStackTrace();
				}
			}
			projectConfig.addSet(classConfig);
		}
		PlainTextConfigurationParser.exportConfiguration(projectConfig, w);
		w.flush();
	}
}