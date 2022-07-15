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

package qspr.metaserver.cs.util;

import java.io.IOException;
import java.io.Serializable;
import java.io.Writer;
import java.util.ArrayList;
import java.util.List;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.metaserver.configurations.LabelWeighting;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration.ScalingType;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.configurations.SupportsInversions;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;

/**
 * Represents a table with labels backed by a DataTable
 * 
 * TODO: Use this class on the client side to create label-tables as well This class an API class, and it should be used everywhere (both clients and server)
 */

public class LabelsTable extends NumericalTable {

	int outputsNumberMultipleProperties;
	List<Integer> optionsNum;
	String implicitConvertedValues [];
	boolean supportsPredicates;
	LabelWeighting labelWeighting; 
	Boolean requireInverse[];

	public LabelsTable() {
	}

	int[] overSample() {
		int oversample[] = new int[getDataSize()];

		return oversample;
	}

	public LabelsTable(DataTable dtValuesOriginal, ModelAbstractConfiguration conf) throws IOException {
		data = dtValuesOriginal;

		supportsPredicates = conf.isSupportPredicates();

		if(conf instanceof MultiLearningAbstractConfiguration)
			labelWeighting = ((MultiLearningAbstractConfiguration)conf).labelWeighting;

		// also values should be provided
		if (data != null && !data.containsColumn(QSPRConstants.VALUE))
			throw new IOException("dtLabels do not contain " + QSPRConstants.VALUE + " column");

		if (conf.getOptions() == null || conf.getOptions().size() == 0) { // to avoid nulls and other uncertainties
			optionsNum = new ArrayList<Integer>();
			optionsNum.add(1);
		} else
			optionsNum = conf.getOptions();

		outputsNumberMultipleProperties = 0;
		for (Integer options : optionsNum) 
			outputsNumberMultipleProperties += options.intValue();

		if (conf.scaleY != null)
			scale = conf.scaleY.get(); // use previous scale
		else 
			if (conf.scaleTypeY != null && conf.scaleTypeY != ScalingType.NONE && !conf.areClassificationData()) 
			{ // we do not have scaling yet -- we create a new one
				if (getDataSize() == 0)
					throw new IOException("There are no data to create scaling!");
				createScaling(conf.scaleTypeY, data.getRowsSize());
				conf.addScaleY(scale);
			}

		if(conf instanceof MultiLearningAbstractConfiguration) { // initialize implicit values substitutions

			MultiLearningAbstractConfiguration config = (MultiLearningAbstractConfiguration) conf;

			if(config.getImplicit() != null) {
				implicitConvertedValues = new String[getColumnSize()];
				data.currentRow = 0;
				//Storing
				Serializable val0 = data.getValue(0);
				Serializable clas0 = data.getValue(1);

				for(int i = 0 ; i < config.getImplicit().length; i++)
					if(config.getImplicit()[i] != null)
					{
						data.setValue(0, config.getImplicit()[i]); // for each missed value
						data.setValue(1, i); // for each class
						String values[] = getScaledValuesString(0); // always the same molecule = 0
						for(int j = 0 ;j < implicitConvertedValues.length; j++) // finding those outputs that are not missed ones
							if(!values[j].equals(LabelsTable.MISSED_VALUES)) {
								String s ="" + config.getImplicit().length+":"; for(Double v : config.getImplicit())s += v + ";";
								if(implicitConvertedValues[j] != null) throw new UserFriendlyException("For implicit: " + s +" value: " + implicitConvertedValues[j]
										+ " is already set instead of " +  values[j] + " length: " + implicitConvertedValues.length);
								implicitConvertedValues[j] = values[j];
							}
					}

				// Restoring
				data.setValue(0, val0);
				data.setValue(1, clas0);
			}

			if(conf instanceof SupportsInversions) {
				for(Integer option: conf.getOptions())
					if(option != 2)
						throw new UserFriendlyException("This method can be used only for data with several properties each classified on two classes only");

				SupportsInversions configuration = (SupportsInversions) conf;

				if(implicitConvertedValues == null)
					throw new UserFriendlyException("This method requires that implicit values will be specified for all classes");

				requireInverse = new Boolean[implicitConvertedValues.length/2];
				for(int i=0; i < requireInverse.length; i++)
					requireInverse[i] = Double.valueOf(implicitConvertedValues[i*2])<0.5;

				outputsNumberMultipleProperties /= 2; // effective number of outputs will be two times smaller
				configuration.setInversions(requireInverse);
			}

		}

	}

	@Override
	public double[] getStd(int size) throws IOException {
		double[] means = getMean(size);
		double std[] = new double[means.length];
		double sizes[] = new double[means.length];

		for (int mol = 0; mol < size; mol++) {
			double v[] = getNonScaledValues(mol);
			for (int i = 0; i < std.length; i++)if(v[i] != MISSED_VALUE) {
				std[i] += (v[i] - means[i]) * (v[i] - means[i]);
				sizes[i]++;
			}
		}

		for (int i = 0; i < std.length; i++)
			std[i] = (float)Math.sqrt(std[i] / sizes[i]);

		return std;
	}

	/**
	 * Return means values of descriptors & value
	 */

	@Override
	public double[] getMean(int size) throws IOException {

		double val[] = getNonScaledValues(0);
		double sizes[] = new double[val.length];

		for(int i=0; i<val.length;i++)
			val[i] = 0;

		for (int mol = 0; mol < size; mol++) {
			double v[] = getNonScaledValues(mol);
			for (int i = 0; i < val.length; i++)
				if(v[i] != MISSED_VALUE) {
					val[i] += v[i];
					sizes[i]++;
				}
		}
		for (int i = 0; i < val.length; i++)
			val[i] /= sizes[i];

		return val;
	}


	@Override
	public double[] getMin(int size) throws IOException {

		double val[] = getNonScaledValues(0);
		for(int i=0; i<val.length;i++)
			val[i] = Float.POSITIVE_INFINITY;

		for (int mol = 0; mol < size; mol++) {
			double v[] = getNonScaledValues(mol);
			for (int i = 0; i < val.length; i++) {
				if (val[i] > v[i] && v[i] != MISSED_VALUE)
					val[i] = v[i];
			}
		}
		return val;
	}

	@Override
	public double[] getMax(int size) throws IOException {

		double val[] = getNonScaledValues(0);
		for(int i=0; i<val.length;i++)
			val[i] = Float.NEGATIVE_INFINITY;

		for (int mol = 0; mol < size; mol++) {
			double v[] = getNonScaledValues(mol);
			for (int i = 0; i < val.length; i++)
				if (val[i] < v[i] && v[i] != MISSED_VALUE)
					val[i] = v[i];
		}
		return val;
	}

	public String getPropertyNameForOutput(int output) {
		return ""+ getPropertyNumberForOutput(output);
	}

	private int getPropertyNumberForOutput(int output) {
		for(int j = 0, num =0 ; j < optionsNum.size(); j++) {
			num += optionsNum.get(j);
			if(output<num)return j;
		}
		return -1;		
	}

	public String getImplicitValue(int molecule, int output) throws IOException {
		if(implicitConvertedValues == null)return null;
		if(requireInverse != null) return getInActiveMap();
		Double implicit[] = (Double [])data.getRow(molecule).getAttachment(QSPRConstants.IMPLICIT_ATTACHMENT);
		if(implicit == null || implicit[getPropertyNumberForOutput(output)] != null) return implicitConvertedValues[output];
		return null; // this value was explicitly requested to be null
	}

	private String getInActiveMap() {
		return "0";
	}

	private String getActiveMap() {
		return "1";
	}


	private String getMap(boolean inverse, boolean on) {
		if(inverse) return on? getInActiveMap() : getActiveMap();
		return on? getActiveMap() : getInActiveMap();
	}


	/**
	 * This function should be used for regression and classification problems with methods which require only one output value
	 */

	public String getOneValueOnlyString(int mol) throws IOException 
	{
		if (optionsNum.get(0) == 1) // Do not scale qualitative labels!
			return scaleValue((Double) data.getValue(mol, 0), 0);
		else
			return ""+ Math.round((Double) data.getValue(mol, 0));
	}

	public int[] countNonMissedNotEmptyClassValues() throws IOException{
		List<String[]> list =new ArrayList<String[]>();
		for (int i=0;i<getDataSize();i++) 
			list.add(getScaledValuesString(i));
		return imbalancedSavedDataCounts(list);
	}

	@Override
	public String[] getScaledValuesString(int molecule) throws IOException {

		String val[] = new String[getColumnSize()];
		double value = (Double) data.getValue(molecule, 0);

		int theValue = (int) Math.rint(value); // the value within the class (property)
		int classID = getClassID(molecule); // we have a value for this class (property)

		int currClass = 0, i = 0;

		for (Integer classOptions : optionsNum) {
			if (classOptions > 1) { // Classification

				if(requireInverse != null && classOptions == 2)  // we have the Taxonomy case here
					val[currClass] = classID == currClass ?  getMap(requireInverse[currClass], theValue == 1): MISSED_VALUES;
				else
					for (int option = 0; option < classOptions; option++) {
						val[i] = classID == currClass ? scaleValue(theValue == option ? 1 : 0, i) : MISSED_VALUES;
						i++;
					}
			}
			else { // Regression
				val[i] = classID == currClass ? getPredicate(molecule) + scaleValue(value, i) : MISSED_VALUES;
				i++;
			}
			currClass++;
		}

		return val;
	}

	/**
	 * Return predicates to be stored in the data file
	 * @param predicate
	 * @return
	 */

	public String getPredicate(int mol) {
		return supportsPredicates? getPredicate((String) data.getRow(mol).getAttachment(QSPRConstants.PREDICATE_ATTACHMENT)):"";
	}


	/**
	 * Maps predicates to those supported with ASNN program
	 * @param predicate
	 * @return
	 */
	String getPredicate(String predicate) {

		if (predicate == null || predicate.equals("=") || predicate.equals("~=") || predicate.equals("~")) 
			return "";

		if (predicate.equals(">") || predicate.equals(">=") || predicate.equals(">>")) 
			return ">";

		if(predicate.equals("<") || predicate.equals("=<") || predicate.equals("<<")) 
			return "<";

		//not mapped to those currently used
		return "";
	}

	/**
	 * Class id of the molecule!
	 * @param molecule
	 * @return
	 */
	public int getClassID(int molecule) {
		return getNumberOfProperties() > 1 ? (int) Math.rint((Double) data.getValue(molecule, QSPRConstants.CLASS)) : 0;
	}

	/**
	 * Counts number of different properties (not outputs!, one property can have many outputs)
	 * 
	 * @return
	 */

	public int getNumberOfProperties() {
		return optionsNum.size();
	}

	/**
	 * Counts number of outputs in the model; one property can have many outputs!
	 * 
	 * @return
	 */

	public int getColumnSize() {
		return outputsNumberMultipleProperties;
	}

	public boolean[] classOutputs() {
		boolean[] out = new boolean[getColumnSize()];
		int i = 0;

		for (Integer classOptions : optionsNum) { 
			for(int j = i; j < i + classOptions; j++) 
				out[j] = classOptions > 1; // this is the classification; allocating correct number of outputs
				i += classOptions;
		}
		return out;
	}

	public int getNumberOfClasserPerProperty(int propertyNumber) throws IOException {
		if (propertyNumber >= optionsNum.size())
			throw new IOException("The requested output does not exist for given configuration");
		return optionsNum.get(propertyNumber);
	}

	public String getMoleculeLabel(int mol) throws IOException {
		int propertyId = getClassID(mol);
		if (getNumberOfClasserPerProperty(propertyId) == 1)
			return propertyId + "_0";
		return propertyId + "_" + getOneValueOnlyString(mol);
	}

	@Override
	public NumericalTable getCopy() {
		LabelsTable t = new LabelsTable();
		t.scale = scale;
		t.optionsNum = optionsNum;
		t.outputsNumberMultipleProperties = outputsNumberMultipleProperties;
		t.implicitConvertedValues = implicitConvertedValues;
		t.requireInverse = requireInverse;
		t.labelWeighting = labelWeighting;
		return t;
	}


	/**
	 * return number of a column that contains maximum double value, i.e. predicted class
	 */

	public static int getPredictedClass(int offset,int classNumber,String [] predictions){
		if(classNumber==1)return 0; // we have regression, i.e. only one class
		return getPredictedClass(offset,classNumber,predictions,-1);
	}

	/**
	 * return number of a column that contains maximum double value, i.e. predicted class
	 * with an exception of a classToExclude 
	 */

	public static int getPredictedClass(int offset,int classNumber,String [] predictions, int classToExclude){

		int clas=-1;
		double val=Double.NEGATIVE_INFINITY;

		for(int i=0;i<classNumber;i++)
			if(i != classToExclude)
			{
				double newval=parseMissedValue(predictions[i+offset]);
				if(val<newval){
					val=newval;
					clas=i;
				}
			}

		return clas;	
	}

	/**
	 *  In some cases, i.e. for intermediate values
	 *  ASNN can provide * instead of real predictions
	 *  In this case we just return 0
	 * @param s
	 * @return
	 */
	public static double parseMissedValue(String s){
		if("*".equals(s))return 0; // in case if there is no real prediction, just return 0
		return Double.parseDouble(s);
	}

	public int[] imbalancedSavedDataCounts(List<String[]> savedValues) {
		int weights[] = new int[getColumnSize()];

		boolean[] outputs = classOutputs(); // indicate whether each particular output is a class or not
		for(String[] vals: savedValues){
			for(int j = 0; j< vals.length; j++)
				if(vals[j] != null && !MISSED_VALUES.equals(vals[j])) {
					if(outputs[j])weights[j] += (Double.valueOf(vals[j]) > 0.5? 1:0 ) ; // for classification we count number of "on" instances
					else
						weights[j]++;
				}
		}

		return weights;
	}

	public double[] getWeights(List<String[]> savedValues) throws IOException{
		if(labelWeighting == null) return null;
		double weights[] = new double [getColumnSize()];

		if(labelWeighting.isGlobalWeighting()) { // calculates frequency and  uses inverse frequency as weight
			int counts[]=imbalancedSavedDataCounts(savedValues);
			for(int i=0;i<weights.length;i++)
				weights[i]=counts[i];

			double sum = 0;
			for(double w:weights)sum += 1./w; sum /= weights.length;
			for(int i=0;i<weights.length;i++)weights[i] = 1./(sum*weights[i]);
			return weights;
		}

		for(int i=0,n=0;i<getNumberOfProperties();i++)
			for(int j=0;j<getNumberOfClasserPerProperty(i); j++) 
				weights[n++] = labelWeighting.getAbsoluteClassWeight(i, j);

		if(labelWeighting.isGlobalNormalization()) {
			double sum = 0;
			for(double w:weights)sum += w; sum /= weights.length;
			for(int i=0;i<weights.length;i++)weights[i] = weights[i]/sum;
		}

		return weights;
	}

	public void saveImplicit(Writer bi) throws IOException {

		if(implicitConvertedValues != null) {
			bi.write("\nimplicit:\t");
			for(String implicit : implicitConvertedValues)
				bi.write(implicit+",");
		}

		if(requireInverse != null) {
			bi.write("\ninverse:\t");
			for(boolean implicit : requireInverse)
				bi.write((implicit ?"1":"0")+",");
		}	
	}

	public boolean isMissedValue(int record, int prop) throws IOException {
		String[]vals= getScaledValuesString(record);
		boolean[] outputs = classOutputs();
		if(vals[prop] == null || MISSED_VALUES.equals(vals[prop]) || (outputs[prop] && Double.valueOf(vals[prop]) < 0.5)) return true; 
		return false;
	}

	public int whichProperty(int record) throws IOException {
		String[]vals= getScaledValuesString(record);
		boolean[] outputs = classOutputs();
		for(int prop=0;prop<vals.length;prop++) {
			if(vals[prop] == null || MISSED_VALUES.equals(vals[prop]) || (outputs[prop] && Double.valueOf(vals[prop]) < 0.5)) continue; 
			return prop;
		}
		throw new IOException("No property has been found for " + record);
	}


}
