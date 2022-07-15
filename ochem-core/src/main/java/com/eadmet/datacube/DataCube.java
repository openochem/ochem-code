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

package com.eadmet.datacube;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

@SuppressWarnings("rawtypes")
public class DataCube<T> {
	
	public List<Subkonto> groupingOrder = new ArrayList<Subkonto>();
	
	public Map<Class, Set> subkontoValues = new LinkedHashMap<Class, Set>();

	private Class<? extends T> clazz;
	private Map<String, T> values = new HashMap<String, T>();
	
	private List<CutDefinition> cutDefinitions = new ArrayList<CutDefinition>();
	
	
	public DataCube(Class<? extends T> clazz) {
		this.clazz = clazz;
	}
	
	public void setupGroupingOrder(Subkonto... subkontos) {
		setupGroupingOrder(Arrays.asList(subkontos));
	}
	
	public void setupGroupingOrder(List<Subkonto> subkontos) {
		this.groupingOrder = subkontos;
		
		cutDefinitions = new ArrayList<CutDefinition>();
		List<Subkonto> subkontoSequence = new ArrayList<Subkonto>();

		cutDefinitions.add(new CutDefinition(subkontoSequence));
		for (int i = 0; i < this.groupingOrder.size(); i++) {
			subkontoSequence.add(subkontos.get(i));
			cutDefinitions.add(new CutDefinition(subkontoSequence));
		}
	}
	
	@SuppressWarnings("unchecked")
	public void addValue(T newValue, Object... addedSubkontoValues) {
			String hash = "";
			String[] cutHashes = new String[cutDefinitions.size()];
			for (int i = 0; i < cutHashes.length; i++)
				cutHashes[i] = "";
			
			for (int i = 0; i < addedSubkontoValues.length; i++)
			{
				Class valueClass = addedSubkontoValues[i].getClass();
				Object value = addedSubkontoValues[i];
				if (valueClass.equals(SubkontoValue.class))
				{
					valueClass = ((SubkontoValue) addedSubkontoValues[i]).subkonto.clazz;
					value = ((SubkontoValue) addedSubkontoValues[i]).value;
				}
				if (!subkontoValues.containsKey(valueClass))
					subkontoValues.put(valueClass, new LinkedHashSet());
				Set entries = subkontoValues.get(valueClass);
				if (!entries.contains(value))
					entries.add(value);
				
				hash += value;
				for (int c = 0; c < cutDefinitions.size(); c++)
					if (cutDefinitions.get(c).containsClass(valueClass))
						cutHashes[c] += value;
			}
		
		if (values.get(hash) == null)
			values.put(hash, newValue);
		else
		{
			T existendValue = values.get(hash);
			existendValue = add(existendValue, newValue);
			values.put(hash, existendValue);
		}
		
		for (int i = 0; i < cutHashes.length; i++) {
			if (cutHashes[i].equals(hash))
				continue;
			T val = values.get(cutHashes[i]);
			
			if (val == null)
				val = getMetricsInstance();
			
			val = add(val, newValue);
			values.put(cutHashes[i], val);
		}
	}
	
	private T add(T destination, T source) {
		if (destination instanceof Double)
			destination = clazz.cast((Double) destination + (Double) source);
		else
			((Metrics) destination).add((Metrics) source);
		return destination;
	}
	
	public Set getSubkontos(String clazz) throws ClassNotFoundException {
		return subkontoValues.get(Class.forName(clazz));
	}
	
	public T getValue(Object... subKontos) {
		String hash = "";
		for (int i = 0; i < subKontos.length; i++)
			hash += subKontos[i];
		T val = values.get(hash);
		return val;
	}
	
	public T getValue(List<Object> subKontos) {
		String hash = "";
		for (int i = 0; i < subKontos.size(); i++)
			hash += subKontos.get(i);
		T val = values.get(hash);
		return val;
	}
	
	public void print() {
		print("", new ArrayList<Object>(), 0);
	}
	
	private void print(String prefix, List<Object> values, int depth) {
		T m = getValue(values);
		if (m != null)
			System.out.println(prefix + "   " + m);
		else
			return;
		
		if (depth >= subkontoValues.size())
			return;
		
		Set currentLevelValues = subkontoValues.get(groupingOrder.get(depth).clazz);
		if (currentLevelValues != null)
		for (Object value : currentLevelValues)
		{
			values.add(value);
			print(prefix + "   ", values, depth + 1);
			values.remove(values.size() - 1);
		}
	}
	
	/**
	 * Get the grouped metrics as a tree-like structure.
	 * Useful for marhalling to XML and gerenating UI reports
	 * @return
	 */
	public DataTree<T> getTree() {
		DataTree<T> root = getTree(new ArrayList<Object>(), 0);
		if (root == null)
		{
			root = new DataTree<T>();
			root.metrics = getMetricsInstance();
			return root;
		}
		root.subkontos = groupingOrder;
		
		return root;
	}
	
	
	private DataTree<T> getTree(List<Object> values, int depth) {
		DataTree<T> node = new DataTree<T>();
		T m = getValue(values);
		if (m != null)
			node.metrics = m;
		else
			return null;
		
		if (!values.isEmpty())
			node.currentGrouping = values.get(values.size() - 1);
		
		if (depth >= groupingOrder.size())
			return node;
		
		Set currentLevelValues = subkontoValues.get(groupingOrder.get(depth).clazz);
		if (currentLevelValues != null)
		for (Object value : currentLevelValues)
		{
			values.add(value);
			DataTree<T> childNode = getTree(values, depth + 1);
			if (childNode != null)
				node.children.add(childNode);
			values.remove(values.size() - 1);
		}
		
		
		
		return node;
	}
	
	private T getMetricsInstance() {
		if (clazz.equals(Double.class))
			return clazz.cast(new Double(0.0));
		try
		{
			return clazz.newInstance();
		} catch (InstantiationException e)
		{
			throw new RuntimeException(e);
		} catch (IllegalAccessException e)
		{
			throw new RuntimeException(e);
		}
	}
}

