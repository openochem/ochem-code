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
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashMap;

import qspr.metaserver.configurations.LinearConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.utils.OSType;

/**
 * 
 * @author itetko
 * Provide specific abstract class layer class for programs developed by IVT
 * 
 */

public abstract class IgorsAbstractServer extends MachineLearningExecutableAbstractServer {

	final static String DATAFILE = "data", CFG = "cfg", OCHEM = "ochem", RESULTS = "results.txt",
			MODEL = "model", OUTPUT = "output", ERROR_PREDICATE = "error:", APPLYFILE = "apply";
	protected String VALUESSEPARATOR = "\t";
	// to put it into prefix

	protected void executeBinary() throws IOException, InterruptedException {

		String[] commands = new String[] { getExeFile(), CFG, DATAFILE };
		executeBinary(commands, OCHEM, 0);
	}

	/**
	 *  Return correct descriptor values (and options converted to descriptors) to be used for any method
	 * @param dtDescriptors
	 * @param molecule
	 * @param configuration
	 * @return
	 * @throws IOException
	 */

	protected double[] getDescriptorValues(DescriptorsTable dtDescriptors, int molecule) throws IOException {
		double vals[]= new double[dtDescriptors.getDescriptorsSize()];
		String descr[]=dtDescriptors.getScaledValuesString(molecule);
		for (int i = 0; i < vals.length; i++)
		{
			descr[i].replace(">", "");
			descr[i].replace("<", "");
			vals[i]=Double.parseDouble(descr[i]);
		}
		return vals;
	}


	/**
	 * Return normalized descriptors
	 * Assumes that table contains them and that they are int/float values
	 * @param coefficient
	 * @param trainSize
	 * @return
	 * @throws IOException 
	 */

	public void addNormaliseCoefficients(DescriptorsTable dtDescriptors,
			LinearConfiguration configuration, int trainSize) throws IOException {

		HashMap<Integer, Double> normalised=new HashMap<Integer, Double> ();
		// Attention -- enumeration in coefficient starts with 1!

		double means[] = dtDescriptors.getMean(trainSize);
		double stds[] = dtDescriptors.getStd(trainSize);

		Double bias=configuration.coefficient.get(0);
		for (Integer num : configuration.coefficient.keySet()){
			if (num > 0){
				bias+=means[num-1]*configuration.coefficient.get(num);
				normalised.put(num, configuration.coefficient.get(num)*stds[num-1]);
			}
		}
		normalised.put(0, bias);

		configuration.normalisedCoefficient=normalised;
	}

	// saves experimental values with or without predicates
	// to be updated to change values for multiple outputs

	protected void saveValues(BufferedWriter writer, LabelsTable dtValues,
			int molecule, int outputs) throws IOException {

		if(dtValues == null || molecule >= dtValues.getDataSize()) // out of the dataset; just add 0s
			writer.append("0");
		else
			writer.append(dtValues.getOneValueOnlyString(molecule));
	}

	public void saveTrainingSetDataforR(String filename, DescriptorsTable dtDescriptors, LabelsTable dtExpValues) throws IOException{

		// create data file
		BufferedWriter writer = getAliasedBufferedWriter(filename);

		// save descriptor names
		for(int i=0;i<dtDescriptors.getDescriptorsSize();i++)
			writer.append("C"+i+VALUESSEPARATOR);
		writer.append("Value\n");

		for(int mol=0;mol<dtExpValues.getDataSize();mol++){

			String v[]=dtDescriptors.getScaledValuesString(mol);
			for(String val:v)
				writer.append(val+VALUESSEPARATOR);
			saveValues(writer, dtExpValues, mol, dtExpValues.getColumnSize());
			writer.append("\n");
		}

		writer.close();
	}

	/**
	 * Saves descriptors and values to the set
	 * If dtExpValues are provided, the number of saved entries is equal to this set
	 * @param filename
	 * @param dtDescriptors
	 * @param dtExpValues
	 * @throws IOException
	 */

	public void saveTrainingSetDataMatrix(String filename, DescriptorsTable dtDescriptors, LabelsTable dtExpValues) throws IOException{

		// create data file
		BufferedWriter writer = getAliasedBufferedWriter(filename);

		int number = dtExpValues == null ? dtDescriptors.getDataSize() : dtExpValues.getDataSize();

		for(int mol=0;mol<number;mol++){

			String v[]=dtDescriptors.getScaledValuesString(mol);
			for(int i=0; i< v.length;i++){
				String val = v[i];
				writer.append(val + (i != v.length - 1 || dtExpValues!=null ? VALUESSEPARATOR : ""));
			}
			if(dtExpValues!=null) saveValues(writer, dtExpValues, mol, dtExpValues.getColumnSize());
			writer.append("\n");
		}

		writer.close();
	}

	/**
	// Saves data from memory to "standard" format supported by ASNN, MLRA, kNN,
	// FSMLR, MLRA programs
	 * 
	 */

	protected void writeSimpleIgorFormatData(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, BufferedWriter writer)
			throws IOException {
		if(dtExpValues == null) throw new IOException("Incorrect funtion is used, for dtExpValues == null indicate the number of outputs values explicitely");
		if(dtDescriptors.getDataSize() < dtExpValues.getDataSize()) throw new IOException("Number of descriptors is less than number of experimental values");
		writeSimpleIgorFormatData(dtDescriptors, dtExpValues, dtExpValues.getColumnSize(), writer,null);
	}


	protected void writeSimpleIgorFormatData(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, BufferedWriter writer, LinkedHashMap<String,Integer> columnsToKeep)
			throws IOException {
		if(dtExpValues == null) throw new IOException("Incorrect function is used, for dtExpValues == null indicate the number of outputs values explicitely");
		if(dtDescriptors.getDataSize() < dtExpValues.getDataSize()) throw new IOException("Number of descriptors is less than number of experimental values");
		writeSimpleIgorFormatData(dtDescriptors, dtExpValues, dtExpValues.getColumnSize(), writer, columnsToKeep);
	}

	protected void writeSimpleIgorFormatData(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, int outputs, BufferedWriter writer, LinkedHashMap<String,Integer> columnsToKeep)
			throws IOException {
		int allMolecules = dtDescriptors.getDataSize();
		int trainingSet = dtExpValues == null? 0: dtExpValues.getDataSize();

		out.println("Number of descriptors: " + (columnsToKeep ==null? dtDescriptors.getDescriptorsSize(): columnsToKeep.size())	+ " all molecules " + allMolecules + " training set: " 
				+ trainingSet);

		if(trainingSet==0)writer.append(allMolecules+ OSType.endLine());
		else
			writer.append(trainingSet
					+ (allMolecules != trainingSet ? " "
							+ (allMolecules - trainingSet) : "") + OSType.endLine());

		boolean cols[]=new boolean[dtDescriptors.getDescriptorsSize()];

		for(int i=0;i<cols.length;i++)
			cols[i]=columnsToKeep==null?true:columnsToKeep.containsKey(dtDescriptors.columnNames.get(i));

		// saving molecules one by one
		for (int mol = 0; mol < allMolecules; mol++){
			writer.append("mol_id" +dtDescriptors.getMoleculeUniqueId(mol)+"_stero"+dtDescriptors.getMoleculeIdStereochem(mol)+VALUESSEPARATOR);

			String vals[]=dtDescriptors.getScaledValuesString(mol);

			for(int i=0;i<vals.length;i++){
				if(cols[i])
					writer.append(vals[i]+VALUESSEPARATOR);
			}

			saveValues(writer, dtExpValues, mol, outputs);
			writer.append(OSType.endLine());
		} 

	}

	/**
	 * Read results in standard format for standard programs, such as KNN, MLRA,
	 * etc. Standard format is mol_xx experimental predicted value
	 * 
	 * @param datafile
	 * @return
	 * @throws IOException
	 */
	protected DataTable readResult(String datafile) throws IOException {

		DataTable dtResult = new DataTable(true);
		dtResult.id = QSPRConstants.PREDICTION_RESULT_COLUMN;
		dtResult.addColumn(QSPRConstants.PREDICTION_RESULT_COLUMN);

		BufferedReader graphics = new BufferedReader(new FileReader(
				getAliasedFileName(datafile)));
		String line;
		while ((line = graphics.readLine()) != null) {
			line = line.trim();
			if (line.startsWith("mol_")) {
				dtResult.addRow();
				String[] prediction = line.split("\\s+");
				dtResult.setValue(QSPRConstants.PREDICTION_RESULT_COLUMN, Double.valueOf(prediction[2]));
			}
		}
		graphics.close();

		return dtResult;

	}


}
