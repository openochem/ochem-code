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

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.apache.commons.lang.StringUtils;


/**
 * Defines the weighting of the labels in the learning process.
 * The weighting can apply to the property as well as to its classes (in case of a classification property)
 * 
 * An exemplary weighting schema:
 * 0.3 LogP
 * 0.7 Ames test
 * 		0.4 active
 * 		0.6 inactive
 * 
 * 
 * @author midnighter
 *
 */
@XmlRootElement
public class LabelWeighting implements Serializable
{
	private static final long serialVersionUID = 1L;

	@XmlElement(name = "property-weight")
	public ArrayList<PropertyWeighting> propertyWeights;

	@XmlElement
	public Boolean useCostMatrix;

	/**
	 * Will normalize parameters inversely proportional number of samples in each class
	 */
	@XmlElement
	public Boolean globalWeighting;

	/**
	 * Normalize everything to sum of weights equals to 1
	 */

	@XmlElement
	public Boolean globalNormalization;

	/**
	 * Create a new Weighting but does not copy the propertyWeights
	 * @param weighting
	 */
	public LabelWeighting(LabelWeighting weighting) {
		useCostMatrix = weighting.useCostMatrix;
		globalNormalization = weighting.globalNormalization;
		globalWeighting = weighting.globalWeighting;
	}

	public LabelWeighting() {
	}

	/**
	 * Get weight of a particular property
	 * @param propertyNumber
	 * @return
	 * @throws IOException 
	 */
	public double getPropertyWeight(int propertyNumber) throws IOException
	{
		if(propertyWeights == null) return 1;
		if(propertyWeights.size() <= propertyNumber || propertyNumber < 0)throw new IOException("No such property " + propertyNumber);
		return propertyWeights.get(propertyNumber).weight;
	}

	public boolean useCostMatrix()
	{
		return useCostMatrix != null && useCostMatrix;
	}

	/**
	 * Make all properties and their classes of equal weight
	 */
	public void equalize()
	{
		if(propertyWeights == null) return;

		for (PropertyWeighting pw : propertyWeights)
		{
			pw.weight = 1;
			if (pw.classesWeights != null)
				for (ClassWeighting cw : pw.classesWeights)
					cw.weight = 1;
		}
	}

	public boolean isGlobalWeighting() {
		return globalWeighting != null && globalWeighting;
	}

	public boolean isGlobalNormalization() {
		return globalNormalization != null && globalNormalization;
	}

	public boolean areWeightsIdentical() {

		double val = propertyWeights.get(0).weight;
		for (PropertyWeighting pw : propertyWeights)
		{
			if(val != pw.weight) return false;
			if (pw.classesWeights == null)continue;
			double valc = pw.classesWeights.get(0).weight;
			for (ClassWeighting cw : pw.classesWeights)
				if(cw.weight != valc) return false;
		}

		return true;
	}

	/**
	 * Do we have more than 1 output? 
	 */
	public boolean hasMultipleOutputs()
	{
		if (propertyWeights == null) return false;

		int outputs = 0;
		for (PropertyWeighting pw : propertyWeights)
			outputs += (pw.classesWeights == null || pw.classesWeights.isEmpty() ? 1 : pw.classesWeights.size());

		return outputs > 1;
	}

	/**
	 * Get relative weight of a particular class in a classification property
	 * @param propertyNumber
	 * @param classNumber
	 * @return
	 * @throws IOException 
	 */
	double getClassWeight(int propertyNumber, int classNumber) throws IOException
	{
		if(propertyWeights == null) return 1;
		if(propertyWeights.size() <= propertyNumber || propertyNumber < 0)
			throw new IOException("No such property " + propertyNumber);
		if(propertyWeights.get(propertyNumber).classesWeights ==null)return 1;
		if (propertyWeights.get(propertyNumber) == null || propertyWeights.get(propertyNumber).classesWeights.size() <= classNumber)
			throw new IOException("No such property weight for property: " + propertyNumber);
		return propertyWeights.get(propertyNumber).classesWeights.get(classNumber).weight;
	}

	/**
	 * Get absolute weight of a particular class in the learning process (also considering the weight of the parent property)
	 * @param propertyNumber
	 * @param classNumber
	 * @return
	 * @throws IOException 
	 */
	public double getAbsoluteClassWeight(int propertyNumber, int classNumber) throws IOException
	{
		if(propertyWeights == null) return 1;
		return getPropertyWeight(propertyNumber) * getClassWeight(propertyNumber, classNumber);
	}

	/**
	 * Add a property or a class with a weight
	 * @param propName the name of the property
	 * @param className the name of the class (can be null)
	 * @param weight the weight in arbitrary scale
	 */
	public void addClass(String propName, String className, double weight)
	{
		if(propertyWeights == null)  propertyWeights =new ArrayList<PropertyWeighting>();

		PropertyWeighting pw = getProperty(propName);
		if (pw == null)
			propertyWeights.add(pw = new PropertyWeighting(propName, weight));
		if (className == null)
			pw.weight = weight;
		else
			if (className != null && pw.getClass(className) == null)
				pw.addClass(className, weight);
	}

	public PropertyWeighting getProperty(String name)
	{
		if(propertyWeights == null)return null;

		for (PropertyWeighting pw : propertyWeights) {
			if (name.equals(pw.name))
				return pw;
		}
		return null;
	}

	public String toString()
	{
		return (propertyWeights == null ? "" : propertyWeights.toString() + (isGlobalNormalization()? " globalNorm": "")) 
				+ (isGlobalWeighting() ? " globalWeighting":"");
	}

	/**
	 * Weighting schema for a particular property
	 */
	public static class PropertyWeighting implements Serializable
	{
		/**
		 * The name of the property. Not required for calculations, but nice to have in the XML
		 */
		@XmlAttribute
		public String name;

		/**
		 * Weights of this property in the multi-property learning
		 */
		@XmlAttribute
		public double weight;

		/**
		 * Weights of of the classes for this property (e.g., "active" can have more weight than "inactive")
		 */
		@XmlElement(name = "property-class")
		public ArrayList<ClassWeighting> classesWeights;

		public PropertyWeighting()
		{

		}

		public PropertyWeighting(String propName, double weight)
		{
			this.name = propName;
			this.weight = weight;
		}

		public ClassWeighting getClass(String name)
		{
			if (classesWeights != null)
				for (ClassWeighting cw : classesWeights)
					if (name.equals(cw.name))
						return cw;
			return null;
		}

		protected void addClass(String className, double weight)
		{
			if (classesWeights == null)
				classesWeights = new ArrayList<LabelWeighting.ClassWeighting>();
			classesWeights.add(new ClassWeighting(className, weight));
		}

		public void setCostMatrix(String class1, String class2, double weight)
		{
			ClassWeighting cw = getClass(class1);
			int ind2 = classesWeights.indexOf(getClass(class2));

			if (cw.costMatrixWeights == null)
				cw.costMatrixWeights = new Double[classesWeights.size()];

			cw.costMatrixWeights[ind2] = weight;
		}

		public double getCostMatrix(String class1, String class2)
		{
			ClassWeighting cw = getClass(class1);
			int ind2 = classesWeights.indexOf(getClass(class2));

			return cw.costMatrixWeights[ind2];
		}

		public String toString()
		{
			String res =  name + " with weight " + weight;
			if (classesWeights != null)
			{
				res += "\n(classes weights: "+classesWeights+")\n";
				if (classesWeights.get(0).costMatrixWeights != null)
				{
					res += "Cost matrix: \n";
					for (int i = 0; i < classesWeights.size(); i++)
						res += classesWeights.get(i).name + ": (" + StringUtils.join(classesWeights.get(i).costMatrixWeights, ", ") + ")\n";
				}
			}

			return res;
		}

		private static final long serialVersionUID = 1L;
	}

	/**
	 * Weighting for a particular class
	 */
	public static class ClassWeighting implements Serializable
	{
		/**
		 * The name of the class (e.g., "active", "inactive"). Not required for calculations, but nice to have in the XML
		 */
		@XmlAttribute
		public String name;

		/**
		 * Weight of of this class
		 */
		@XmlAttribute
		public double weight;

		@XmlElement
		public Double[] costMatrixWeights;

		public ClassWeighting()
		{

		}

		public ClassWeighting(String className, double weight)
		{
			this.name = className;
			this.weight = weight;
		}

		public ClassWeighting(String className, double weight, Double[] costMatrixRow)
		{
			this.name = className;
			this.weight = weight;
			this.costMatrixWeights = costMatrixRow;
		}

		public String toString()
		{
			return name + "*" + weight;
		}

		private static final long serialVersionUID = 1L;
	}


	public void addNewProperty(PropertyWeighting pwNew) {
		if(propertyWeights == null)  propertyWeights =new ArrayList<PropertyWeighting>();
		propertyWeights.add(pwNew);		
	}

}

