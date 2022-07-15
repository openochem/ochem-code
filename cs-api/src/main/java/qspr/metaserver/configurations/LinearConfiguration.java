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

package qspr.metaserver.configurations;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import javax.xml.bind.annotation.XmlRootElement;

import com.eadmet.utils.NumericalValueStandardizer;


@XmlRootElement(name = "linear-configuration")
abstract public class LinearConfiguration extends ModelAbstractConfiguration {

	private final static String BIAS = "Bias";
	private final static String ERROR = "Error";
	private final static String DESCRIPTOR = "Descriptor";
	private static final long serialVersionUID = 3L;

	/**
	 * The first coefficient is the bias and should be always present
	 * All other coefficients start with 1, i.e. 
	 * coefficient[1,VAL] contains the coefficient for the first variable (not for the 0!)
	 */

	public Map<Integer, Double> coefficient = new HashMap<Integer, Double>();
	public Map<Integer, Double> normalisedCoefficient = new HashMap<Integer, Double>();

	public Double minVal,maxVal;
	public Boolean limitRange;

	@Override
	public boolean isTrainingConfiguration()
	{
		return (coefficient==null || coefficient.isEmpty()) && savedmodel == null;
	}

	public boolean isModelSaved() {
		return !isTrainingConfiguration();
	}
	
	@Override
	public List<Integer> getSelectedDescriptors()
	{
		List<Integer> descNums = new ArrayList<Integer>();
		for (Integer num : coefficient.keySet())
			if (num != 0) // coefficient 0 corresponds to the intercept
				descNums.add(num - 1);

		Collections.sort(descNums);
		return descNums;
	}

	public String theEquation(List<String> descriptors, boolean normalized){
		StringBuffer equation = new StringBuffer();

		if(coefficient==null || coefficient.size()==0)return null;

		Map<String, Serializable> eq=getEquation(descriptors, normalized);

		for(String descName: eq.keySet()){

			if(descName.equals(DESCRIPTOR))continue;
			if(descName.equals(ERROR))return (String)eq.get(descName);

			double value=(Double)eq.get(descName);

			if(descName.equals(BIAS))equation.append(" Y = " + NumericalValueStandardizer.getSignificantDigits(value));
			else{
				equation.append((value>0?" + ":" - ") + NumericalValueStandardizer.getSignificantDigits(value>0?value:-value) 
						+ "*" + descName);
			}
		}

		if(normalized)equation.append("\n("+(coefficient.size()-1)+" variables in equation)");

		return equation.toString();
	}

	public Map<String, Serializable> getEquation(List<String> descriptors, boolean normalized)
	{
		Map<String, Serializable> equation = new LinkedHashMap<String, Serializable>();
		equation.put("Descriptor", "Coefficient");

		Map<Integer, Double> coeff= normalized?normalisedCoefficient:coefficient;
		if(coeff==null || coeff.size()==0)
			equation.put("equation", "Model was not saved");
		else
		{
			Map<Integer, Double> sorted = normalized ? sortHashMapAbsoluteValues(normalisedCoefficient) : coefficient;

			equation.put(BIAS, Double.parseDouble(NumericalValueStandardizer.getSignificantDigits(coeff.get(0))));

			try
			{
				for (Integer descNum : sorted.keySet())
					if (descNum != 0){
						double val=Double.parseDouble(NumericalValueStandardizer.getSignificantDigits(coeff.get(descNum)));
						if(val!=0)equation.put(descriptors.get(descNum - 1),val);
					}
			}
			catch (IndexOutOfBoundsException e)
			{
				equation = new LinkedHashMap<String, Serializable>();
				equation.put(ERROR, "Invalid model, it refers to unexistent descriptors: coeff:"+coeff.size()
						+" desc: "+descriptors.size()+ equation+e.getLocalizedMessage());
			}
		}

		return equation;
	}

	/*
	 * Provides sorting of values in the decreasing order to identify the most significant variables
	 */
	private Map<Integer, Double> sortHashMapAbsoluteValues(Map<Integer, Double> input)
	{
		Map<Integer, Double> tempMap = new HashMap<Integer, Double>();
		for (Integer wsState : input.keySet()){
			tempMap.put(wsState,Math.abs(input.get(wsState))); // by absolute values
		}

		List<Integer> mapKeys = new ArrayList<Integer>(tempMap.keySet());
		List<Double> mapValues = new ArrayList<Double>(tempMap.values());
		HashMap<Integer, Double> sortedMap = new LinkedHashMap<Integer, Double>();
		TreeSet<Double> sortedSet = new TreeSet<Double>(mapValues);
		Object[] sortedArray = sortedSet.toArray();
		int size = sortedArray.length;
		for (int i=size-1; i>=0; i--){ // we sort in the decreasing order of absolute values
			sortedMap.put(mapKeys.get(mapValues.indexOf(sortedArray[i])), 
					(Double)sortedArray[i]);
		}
		return sortedMap;
	}

	public String toString(){
		return minVal != null ? "minValue="+NumericalValueStandardizer.getSignificantDigitsStr(minVal,3)+" maxValue="+
				NumericalValueStandardizer.getSignificantDigitsStr(maxVal,3):" " + super.toString();
	}

	@Override
	public boolean isSupportRegression() {
		return true;
	}
}



