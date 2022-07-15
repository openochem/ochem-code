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

package qspr.metaserver.cs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.Writer;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.utils.NumericalValueStandardizer;

import qspr.metaserver.configurations.EarlyStopping;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.configurations.SupportsInversions;
import qspr.metaserver.configurations.SupportsOneOutputOnly;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

public abstract class MultiLearningAbstractServer extends MachineLearningExecutableAbstractServer{

	final static String WEIGHTS = "weights.txt";
	private final static String IDS = "molecules.ids.csv";
	final static String GRAPH = "graph.txt";
	protected String DELIMITER = ",";

	/**
	 * 
	 * @param filename
	 * @param dtDescriptors 
	 * @param dtExpValues - can be null; saving only descriptors
	 * @param configuration
	 * @return - number of correctly saved records
	 * @throws Exception
	 */
	protected int saveAggregatedData(String filename, DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration configuration) throws Exception {
		// saving all values
		BufferedWriter bw = getAliasedBufferedWriter(filename);
		int savedRecords = saveAggregatedDescriptorsAndValues(dtDescriptors, dtExpValues, bw, configuration);
		bw.close();

		LineNumberReader read = getAliasedLineNumberReader(filename);
		String s = null;
		while((s = read.readLine()) != null) 
			if(s.contains("null")) {
				read.close();
				throw new IOException("null values were saved in the file. Calculations with "+ configuration.getInformativeName() + " is not possible.");
			}
		read.close();

		return savedRecords;
	}

	@Override
	public boolean isCritical(String message) {
		if(message.contains("constant value"))
			throw new CriticalException("For (one of) target property all values were identical in a validation fold. This is not allowed. Increase dataset size or decrease number of folds.");
		if(message.contains("only missed values are in column"))
			throw new CriticalException(message+ " - check that your data have meanigful values to predict");
		if(message.contains("n range of some loaded parameters"))
			throw new CriticalException(message+ " - check that your conditions are covered by the training set");
		return false;
	}

	@Override
	final protected DataTable applyModelBatch(DescriptorsTable dtDescriptors, ModelAbstractConfiguration receivedConf) throws Exception
	{
		MultiLearningAbstractConfiguration conf = (MultiLearningAbstractConfiguration) receivedConf;
		DataTable tab =  applyModelMulti(dtDescriptors,conf);

		if(conf.getImplicit() != null)
			for(int i=0;i<tab.getRowsSize();i++)
				if(tab.getRow(i).getAttachment(QSPRConstants.IMPLICIT_ATTACHMENT) == null)
					tab.getRow(i).addAttachment(QSPRConstants.IMPLICIT_ATTACHMENT, conf.getImplicit()); // we add for all new predictions too

		return tab;
	}

	abstract protected DataTable applyModelMulti(DescriptorsTable dtDescriptors, MultiLearningAbstractConfiguration receivedConf) throws Exception;
	abstract protected void saveOneRow(int mol, DescriptorsTable dtDescriptors, String[] values, Writer bw) throws IOException;

	final public int saveAggregatedDescriptorsAndValues(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, Writer bw, ModelAbstractConfiguration configration) throws IOException {

		Writer bi = null;

		if(dtExpValues == null) { 
			if (bw != null)bi = getAliasedBufferedWriter(IDS+".apply");
			int n = saveDescriptors(dtDescriptors, bw, bi); 
			if (bi != null) bi.close();
			return n;
		}

		Map<String,ArrayList<Map.Entry<Integer,String>>> val = new LinkedHashMap<String,ArrayList<Map.Entry<Integer,String>>> ();
		ArrayList<Map.Entry<Integer,String>> vals ;
		int outputs = dtExpValues.getColumnSize();

		if(configration instanceof SupportsOneOutputOnly && outputs ==2 && configration.areClassificationData())
			outputs = 1;

		List<String[]> savedValues = new ArrayList<String[]>();

		int count = 0, molecules = dtExpValues.getDataSize();

		try {
			if(bw != null)bi = getAliasedBufferedWriter(IDS);

			// collection of data according to hash values by descriptors and conditions
			for (int i = 0; i < molecules; i++)
				if(!dtDescriptors.getRawRow(i).isError())
				{
					String hashVal = dtDescriptors.getMoleculeMD5(i, true);
					if(val.containsKey(hashVal)) // val contains for all molecules
						vals = val.get(hashVal);  // values for this molecule, i.e. output 10:0.1 - for output 10 we have value 0.1 
					else
					{
						vals = new ArrayList<Map.Entry<Integer,String>>();
						val.put(hashVal,vals);
					}

					String values[] = dtExpValues.getScaledValuesString(i);
					for(int n = 0; n < outputs; n++)
						if(!values[n].equals(LabelsTable.MISSED_VALUES))
							vals.add(new AbstractMap.SimpleEntry<Integer,String>(n,values[n])); // all non missed values are stored for each descriptor hash
				}
				else
					throw new IOException("No errors are allowed in this function!");

			setStatus("Unique molecules: " + val.size() + " out of total: " + molecules);


			Integer array[] = new Integer[molecules];
			for (int i = 0; i < array.length; i++) array[i] = i;
			int seed = configration.getSeed();
			Collections.shuffle(Arrays.asList(array), new Random(seed));
			for(int i = 0; i< array.length;i++)if(array[i] == 0) { // required to print headers 
				array[i] = array[0];
				array[0] = 0;
				break;
			}

			String  strOriginal[] = new String[outputs];

			Set<String> earlyStopping = new HashSet<String>();

			if(configration instanceof EarlyStopping) {
				EarlyStopping conf = (EarlyStopping) configration;
				if(!conf.shuffleAllData()) {
					double fraction = conf.getEarlyStoppingFraction();
					Set<String> ha = val.keySet();
					Object hashes[]= ha.toArray();
					Collections.shuffle(Arrays.asList(hashes), new Random(seed+1));
					for(int i = 0; i<hashes.length*fraction && i<hashes.length;i++)
						earlyStopping.add((String)hashes[i]);
					setStatus("Using "+earlyStopping.size() + " mols out of " + hashes.length + " for early stopping" );
				}
			}

			for(Integer i: array) { // selecting next molecule for training subset
				String hashVal = dtDescriptors.getMoleculeMD5(i, true); // selection by descriptors and conditions
				if(earlyStopping.contains(hashVal))continue;
				count += saveMol(i, dtDescriptors, val, outputs, strOriginal, dtExpValues, savedValues, bw, bi);
			}

			if(earlyStopping.size() > 0)
				for(Integer i: array) { // selecting next molecule for early stopping subset
					String hashVal = dtDescriptors.getMoleculeMD5(i, true); // selection by descriptors and conditions
					if(!earlyStopping.contains(hashVal))continue;
					count += saveMol(i, dtDescriptors, val, outputs, strOriginal, dtExpValues, savedValues, bw, bi);
				}

		}catch(Exception e) {
			throw e;
		}
		finally {
			if(bi != null) {
				dtExpValues.saveImplicit(bi);
				bi.write("\n");
				bi.close();
			}
		}

		saveWeights(dtExpValues.getWeights(savedValues));

		return count;
	}

	void saveAnnotation(int i, DescriptorsTable dtDescriptors, String  strVals[], String  strOriginal[], Writer bi) throws IOException {
		String hashVal = dtDescriptors.getMoleculeMD5(i, true); // selection by descriptors and conditions

		bi.write(""+ i + "," + 
				dtDescriptors.getRawRow(i).attachments.get(QSPRConstants.SMILES_ATTACHMENT ) + ","  + 
				dtDescriptors.getRawRow(i).attachments.get(QSPRConstants.MOLECULE_ID_STEREOCHEM ) + "," +
				(dtDescriptors.getRawRow(i).attachments.get(QSPRConstants.EXTERNAL_ID) != null ? 
						dtDescriptors.getRawRow(i).attachments.get(QSPRConstants.EXTERNAL_ID) + "," : "")+
				dtDescriptors.getRawRow(i).attachments.get(QSPRConstants.INCHIKEYS)+"," +
				hashVal + ","
				);
		if(strVals !=null) {
			DescriptorsTable.saveOneRowValuesSVM(strVals,bi); 
			bi.write(",/,");
			DescriptorsTable.saveOneRowValuesSVM(strOriginal,bi);
		}
		bi.write("\n");
	}

	int saveMol(int i, DescriptorsTable dtDescriptors, Map<String,ArrayList<Map.Entry<Integer,String>>> val,
			int outputs, String  strOriginal[], LabelsTable dtExpValues, List<String[]> savedValues, Writer bw, Writer bi) throws IOException {

		String hashVal = dtDescriptors.getMoleculeMD5(i, true); // selection by descriptors and conditions
		if(!val.containsKey(hashVal))return 0; // already was saved...

		String  strVals[] = new String[outputs]; // original is required each time to properly store numbers of records

		ArrayList<Map.Entry<Integer,String>> vals  = val.get(hashVal);
		for(int n = 0; n < outputs; n++){
			Map.Entry<Integer,String> foundValue = null;
			for(Map.Entry<Integer,String> p:vals){
				if(p.getKey() != n)continue; // we need values for this output
				foundValue = p;
				break;
			}

			if(foundValue != null){
				strOriginal[n] = strVals[n] = foundValue.getValue();
				vals.remove(foundValue);
			}
			else { //
				strVals[n] = dtExpValues.getImplicitValue(i, n);
				strOriginal[n] = null;
			}
		}
		if(vals.size() == 0)val.remove(hashVal); // no new values to save for this descriptor
		if(bw != null) {
			if(bi != null)saveAnnotation(i, dtDescriptors, strVals, strOriginal, bi); 
			saveOneRow(i, dtDescriptors, strVals, bw);
			savedValues.add(strVals);
		}
		return 1;
	}

	private void saveWeights(double val[]) throws IOException{
		if(val == null) return;
		BufferedWriter writer  = getAliasedBufferedWriter(WEIGHTS);
		for(int i = 0 ; i< val.length; i++)
			writer.write( (i==0?"":",") + NumericalValueStandardizer.getSignificantDigitsDouble(val[i], NumericalValueStandardizer.SIGNIFICANT_DIGITS));
		writer.write("\n");
		writer.close();
	}

	private int saveDescriptors(DescriptorsTable dtDescriptors, Writer bw, Writer bi) throws IOException {
		if(bw != null)
			for (int i = 0; i < dtDescriptors.getDataSize(); i++) {
				if(!dtDescriptors.getRawRow(i).isError())
					saveOneRow(i,dtDescriptors, null, bw);
				saveAnnotation(i, dtDescriptors, null, null, bi);
			}
		return dtDescriptors.getRawData().getRowsNoErrorsSize();
	}

	public DataTable readResultValues(String filename, ModelAbstractConfiguration conf) throws IOException{
		return readResultValues(filename, conf, 0);
	}

	public DataTable readResultValues(String filename, ModelAbstractConfiguration conf, int startPosition) throws IOException{
		List<Integer> optionsNumber = conf.optionsNumber;

		BufferedReader br = exeRunner.getAliasedBufferedReader(filename);

		DataTable dtResult = new DataTable(true);
		dtResult.id = QSPRConstants.PREDICTION_RESULT_COLUMN;

		int optionsSize = optionsNumber == null ? 1: optionsNumber.size();

		Boolean inverse[] = null;

		if(conf instanceof SupportsInversions) {
			SupportsInversions config = (SupportsInversions) conf;
			inverse = config.getInversions(); // do we have them?!
			if(inverse != null)
				setStatus("Read results using inversions");
			else
				setStatus("Read results");
		}

		int std = 0;
		String line;
		while ((line = br.readLine()) != null)
			try{
				
				line = line.trim();
				
				String[] predictions = line.split(DELIMITER);

				if(line.contains(QSPRConstants.PREDICTION_RESULT_COLUMN) && line.contains(QSPRConstants.DM)) {
					std = (predictions.length - startPosition)/2; 
				}
				if(line.contains(QSPRConstants.PREDICTION_RESULT_COLUMN) || line.contains("\"x\"")) continue; // Second is for R---

				dtResult.addRow();
				if(line.toUpperCase().contains(QSPRConstants.ERROR)) {
					dtResult.getCurrentRow().setError("prediction failed");
					continue;
				}

				// TODO to refactor, only offset should be used
				for (int output = 0, offset = startPosition, classNumber; output < optionsSize; output++, offset += classNumber){
					String name= optionsSize == 1 ? "" : "" + output;

					classNumber= optionsNumber == null ? 1 :  optionsNumber.get(output);

					if(classNumber > 2) {
						int predictedClass = LabelsTable.getPredictedClass(offset,classNumber,predictions);
						dtResult.setValue(QSPRConstants.PREDICTION_RESULT_COLUMN + name, predictedClass); // class for multi-class classification
						if(std>0)addSTD(dtResult, name, predictions[offset + std]);
						continue;
					}

					if(classNumber>1 & inverse != null){ //we have a classification task on two classes
						Double value = Double.valueOf(predictions[startPosition + output]);
						value = correct(inverse[output] ? 1.-value : value);
						dtResult.setValue(QSPRConstants.PREDICTION_RESULT_COLUMN + name, value);
						if(std>0)addSTD(dtResult, name, predictions[offset + std]);
						continue;
					}

					if(classNumber>1)
						dtResult.setValue(QSPRConstants.PREDICTION_RESULT_COLUMN + name, correct(1.-Double.parseDouble(predictions[offset])));
					else
						dtResult.setValue(QSPRConstants.PREDICTION_RESULT_COLUMN + name, Double.parseDouble(predictions[offset]));

					if(std>0)addSTD(dtResult, name, predictions[offset + std]);

				}
			}catch(Throwable e) {
				dtResult.getCurrentRow().setError(e.getMessage());
			}

		br.close();
		return dtResult;
	}

	private void addSTD(DataTable dtResult, String name, String val) {
		//dtResult.setValue(QSPRConstants.DM +  name + ":STD", val);
	}

}
