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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.NoDescriptors;
import qspr.metaserver.configurations.ValidationConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration.ScalingType;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.CompactDataRow;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.utils.NumericalValueStandardizer;
import com.eadmet.utils.OCHEMUtils;

/**
 * Represents a table with molecular descriptors backed by a DataTable.
 * Provides correct data normalization
 * 
 */

public class DescriptorsTable extends NumericalTable {

	static final float QUALITATIVE_OPTION_THRESHOLD = 0.000001f;

	Map<Integer, Integer> qualitativeDescriptors;

	public List<String> columnNames;

	public boolean considerSmiles = false;

	/**
	 * Counts only number of columns that are descriptors and are not conditions
	 */

	int onlyDescriptors = 0; 

	private DescriptorsTable() {
	}

	public int getConditionsSize() {
		return getColumnSize() - onlyDescriptors;
	}

	public DescriptorsTable(DataTable dtDescriptors, ModelAbstractConfiguration configuration, int trainingSetSize) throws IOException {

		considerSmiles = configuration instanceof NoDescriptors
				|| (configuration instanceof ValidationConfiguration && ((ValidationConfiguration)configuration).taskConfiguration instanceof NoDescriptors)
				; // individual SMILES 

		// we are in the apply model and should use previous cfg
		if (trainingSetSize == 0) {
			loadTable(dtDescriptors, configuration.qualitativeDescriptors);
			if (configuration.scaleX != null)
				scale = configuration.scaleX.get();
		} else {
			// we have to create it and save it
			loadTable(dtDescriptors, calculateQualitativeDescriptorsMap(dtDescriptors));
			if (configuration != null) {
				configuration.qualitativeDescriptors = getQualitativeDescriptorsMap(); // store
				if (configuration.scaleTypeX != null && configuration.scaleTypeX != ScalingType.NONE) {
					createScaling(configuration.scaleTypeX, trainingSetSize);
					configuration.addScaleX(scale);
				}
			}
		}

		for (onlyDescriptors = 0; onlyDescriptors < dtDescriptors.getColumnsSize(); onlyDescriptors++) 
			if(data.getColumnAttachment(data.getColumn(onlyDescriptors), QSPRConstants.IS_CONDITION_COLUMN) != null)break;

	}

	@Override
	public int getColumnSize() {
		return  columnNames.size();
	}

	private void loadTable(DataTable dtDescriptorsOriginal, Map<Integer, Integer> qualitativeDesc) throws IOException {

		data = dtDescriptorsOriginal;

		// if we have qualitative descriptors, we need to create additional mappings
		if (qualitativeDesc != null) {
			qualitativeDescriptors = qualitativeDesc;
			columnNames = new ArrayList<String>();
			for (int i = 0; i < data.getColumnsSize(); i++) {
				if (qualitativeDescriptors.containsKey(i))
					for (int j = 0; j < qualitativeDescriptors.get(i); j++) 
						columnNames.add(data.getColumn(i) + "_" + j);					
				else 
					columnNames.add(data.getColumn(i));
			}
		} else {
			qualitativeDescriptors = new HashMap<Integer, Integer>();
			columnNames = data.getColumns();
			if(columnNames == null) columnNames = new ArrayList<String>();
		}
	}

	public double getRawColumnValue(int mol, int option) {
		return (Double) data.getValue(mol, option);
	}

	public String getRawColumnName(int option) {
		return data.getColumn(option);
	}

	public int getRawColumnsNumber() {
		return data.getColumnsSize();
	}

	public int optionSize(int option) {
		return qualitativeDescriptors.get(option) == null ? 1 : qualitativeDescriptors.get(option);
	}

	public String getDescriptorName(int descriptor) {
		return columnNames.get(descriptor);
	}

	public int getDescriptorsSize() {
		return getColumnSize();
	}

	public int getNonZerosDescriptorsSize(int rows) throws IOException {
		int maxv = 0, j = 0;
		for(int i=0;i<rows;i++){
			double values[] = getScaledValues(i);
			for(j = values.length - 1; j > maxv && values[j] == 0 ; j--);
			if(j > maxv)maxv = j;
		}
		return maxv + 1;
	}

	static public void saveOneRowValuesSVM(String[] values,  Writer out) throws IOException {
		if(values != null)
			for(int i= 0; i < values.length; i++) {
				out.write(values[i] != null?values[i]:"");
				if( i != values.length -1 )out.write(",");			
			}else
				out.write("0");
	}

	public void saveOneRowSVM(int mol, String[] values,  Writer out) throws IOException {

		saveOneRowValuesSVM(values, out);

		double desc[] = getScaledValues(mol);
		for (int i=0; i< desc.length; i++)
		{
			double v = desc[i];
			if (v != 0)
				out.write(" " + (i+1) + ":" +  NumericalValueStandardizer.getSignificantDigits(v));
		}

		if( mol == 0 && desc[desc.length - 1] == 0) // we need to add zero just to indicate length of data for SVM and fix bug
			out.write(" " + desc.length + ":0");

		out.write("\n");
	}

	public String getMoleculeMD5(int mol, boolean useAlsoConditions) throws IOException{

		if(!data.isCompactRowFormat())throw new IOException("Descriptors data are not in the compact format.");

		String val[] = getScaledValuesString(mol);
		StringBuffer b = new StringBuffer(getDescriptorsSize() * 10);
		for (int i = 0; i < val.length && (useAlsoConditions || i < onlyDescriptors); i++) 
			b.append(val[i]+"_");

		if(considerSmiles) { // special case for data with SMILES and also augmentations
			String inchies = (String)data.getRow(mol).getAttachment(QSPRConstants.INCHIKEYS);
			if(inchies == null)
				throw new IOException("No Inchies were provided for mol: " + mol+ " RECORDID: " + data.getRow(mol).getAttachment(QSPRConstants.RECORD_ID_ATTACHMENT));
			b.append(inchies);
		}

		//return b.toString();  //TODO debug!!!
		return OCHEMUtils.getMD5(b.toString());	
	}

	@Override
	public String[] getScaledValuesString(int molecule) throws IOException{

		if (molecule >= getDataSize())
			throw new IOException("molecule " + molecule + " > datasize" + getDataSize());

		float [] values=((CompactDataRow)data.getRow(molecule)).toArray();  // just raw values 

		String val[] = new String[columnNames.size()];

		for (int i = 0, n=0; n < val.length; i++) { // values length can be smaller than that of number of descriptors (in case if we have qualitative conditions)

			float value = i < values.length ? values[i] : 0;

			if (!qualitativeDescriptors.containsKey(i)) {
				val[n] = "" + floatScaleValue(value, n);
				n++;
			} else {
				int option = (int) Math.rint(value);
				int maxOpt = qualitativeDescriptors.get(i);
				if (Math.abs(value - option) > QUALITATIVE_OPTION_THRESHOLD || option >= maxOpt)
					throw new IOException("Wrong opt value " + value + " maxOpt=" + maxOpt);

				for (int j = 0; j < maxOpt; j++) {
					val[n] = "" + floatScaleValue(option == j ? 1 : 0, n);
					n++;
				}
			}
		}

		return val;
	}

	public Integer getMoleculeUniqueId(int mol) {
		return (Integer) data.getRow(mol).getAttachment(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM);
	}

	public Integer getMoleculeIdStereochem(int mol) {
		return (Integer) data.getRow(mol).getAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM);
	}

	private Map<Integer, Integer> calculateQualitativeDescriptorsMap(DataTable dtDescriptors) {
		Map<Integer, Integer> qualDescriptors = new HashMap<Integer, Integer>();
		System.out.println("Filling in the qualitative descriptors");

		for (int col = 0; col < dtDescriptors.getColumnsSize(); col++)
			if (dtDescriptors.getColumnAttachment(dtDescriptors.getColumn(col), QSPRConstants.IS_QUALITATIVE_COLUMN) != null) {
				List<Serializable> optionList = dtDescriptors.columnToListUnique(col);
				int maxOption = 0;
				for (Serializable option : optionList) {
					double value = (Double) option;
					if (Math.rint(value) >= maxOption)
						maxOption = (int) Math.rint(value)+1;
				}
				//				if(maxOption<=1) // this is error; user may not be aware that there is  no variation in this descriptor
				//					throw new UserFriendlyException(
				//							"There is only one option for qualitative condition \""+dtDescriptors.getColumn(col)+"\"\n"+
				//							"Eliminate this condition from the list of descriptors.");
				System.out.println("Number of options for " + dtDescriptors.getColumn(col) + ": " + maxOption);

				qualDescriptors.put(col, maxOption);
			}

		System.out.println("" + qualDescriptors.size() + " descriptors are qualitative");

		return qualDescriptors.size() == 0 ? null : qualDescriptors;
	}

	private Map<Integer, Integer> getQualitativeDescriptorsMap() {
		return qualitativeDescriptors;
	}


	public DescriptorsTable getSliceDescriptorsTable(int fromIndex, int toIndex) {
		DescriptorsTable table = (DescriptorsTable) getCopy();
		table.data = data.getSlice(fromIndex, toIndex);
		return table;
	}

	@Override
	public DescriptorsTable getCopy() {
		DescriptorsTable t = new DescriptorsTable();
		t.data = data;
		t.qualitativeDescriptors = qualitativeDescriptors;
		t.columnNames = columnNames;
		t.scale = scale;
		t.considerSmiles = considerSmiles;
		t.onlyDescriptors = onlyDescriptors;
		return t;
	}

	public DescriptorsTable compressDescriptorsForApply(int[] compressing) throws IOException{
		DescriptorsTable res = getCopy();
		res.data = data.getEmptyCopy();
		HashMap<String,Integer> hash = new HashMap<String,Integer>();
		int n =0;
		for(int i=0;i<getDataSize();i++) {
			String newhash=getMoleculeMD5(i,true);
			if(!hash.containsKey(newhash)) {
				res.data.addRow(data.getRow(i));
				hash.put(newhash, n++);
			}
			compressing[i] = hash.get(newhash);
		}
		return res;
	}

	public Map<String, Integer[]> groupMoleculesByMD5Hashes(LabelsTable values) throws IOException {
		Map<String,Integer[]> hashes = new HashMap<String,Integer[]>(); 
		int size = values.getNumberOfProperties();
		for (int i = 0; i < values.getDataSize(); i++){
			String hash = getMoleculeMD5(i, true);
			if(!hashes.containsKey(hash)) { 
				Integer val[] = new Integer[size];
				for(int j=0;j<size;j++)val[j]=0;
				hashes.put(hash,val);
			}
			hashes.get(hash)[values.getClassID(i)]++;
		}
		return hashes;
	}

	public int getUniqueMoleculesNumber() {
		HashSet<Integer> r = new HashSet<Integer>();
		data.reset();
		while(data.nextRow()){
			r.add((Integer) data.getCurrentRow().getAttachment(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM));
		}

		return r.size();
	}

}
