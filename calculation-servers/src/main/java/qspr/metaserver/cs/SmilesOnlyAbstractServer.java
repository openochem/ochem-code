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
import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ThreadLocalRandom;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.OSType;

import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.configurations.NoDescriptors;
import qspr.metaserver.configurations.Stereochemistry;
import qspr.metaserver.configurations.SupportsOneOutputOnly;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.protocol.Task.DebugLevel;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

abstract public class SmilesOnlyAbstractServer extends MultiLearningAbstractServer{

	protected String executablePython = null;
	final private String output = "output";

	final String CFG = "config.cfg";
	final String DATAFILE = "train.csv";
	final String APPLYFILE = "apply.csv";
	String MODEL = "model.tar";
	final String PREDICTIONS = "results.csv";

	int startPositionResults = 1;

	abstract void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean traningMode, boolean forceCPU) throws IOException, InterruptedException;

	private String getFilter() {
		File f = new File(getExeFile());
		if(!f.isDirectory())
			f = f.getParentFile();
		String s = f.getParentFile().getParent();
		s += "/tools/utils/filter.py";
		return s;
	}

	@Override
	protected void saveOneRow(int mol, DescriptorsTable dtDescriptors, String[] values,  Writer writer) throws IOException {
		String smile = (String) dtDescriptors.getRawData().getRow(mol).getAttachment(QSPRConstants.SMILES_ATTACHMENT);
		writer.write(smile);
		if(dtDescriptors.getDescriptorsSize()>0) {
			String []desc = dtDescriptors.getScaledValuesString(mol);
			for(int i=0; i <desc.length; i++)
				writer.append(","+ (desc[i] == null ? "": desc[i]));
		}

		if(values != null)
			for(int i=0; i <values.length; i++)
				writer.append(","+ (values[i] == null ? "": values[i]));
		writer.append(OSType.endLine());
	}

	/**
	 * Provides standard conversion of SMILES or just raw SMILES (or data)
	 * @param mol
	 * @param dtDescriptors
	 * @param convertBeforeSaving
	 * @return
	 * @throws IOException
	 */


	String getQualifiedName( ModelAbstractConfiguration receivedConf){
		return (supportedTaskType.equals(receivedConf.getInformativeName()) ?"":" " + receivedConf.getInformativeName());
	}

	protected CONDA getCondaEnvironment(){
		return CONDA.RDKIT;
	}

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws Exception {

		int initialSize = dtDescriptors.getDataSize();

		System.out.println("mols="+dtDescriptors.getDataSize() + " noerr="+dtDescriptors.getRawData().getRowsNoErrorsSize());

		setStatus("Starting " + getQualifiedName(receivedConf));
		saveConfig(receivedConf, dtDescriptors, true, false); // order is important (before descriptors) to initialize DeepChem method

		saveAggregatedData(getAliasedFileName(DATAFILE), dtDescriptors, dtExpValues, receivedConf);
		if(debug > DebugLevel.NONE)
			saveAggregatedData(getAliasedFileName(APPLYFILE), dtDescriptors, null, receivedConf);

		setStatus("Running " + getQualifiedName(receivedConf));

		executablePython = runPythonWithConda(getCommands(receivedConf,true), MODEL, getCondaEnvironment());

		receivedConf.storeModel(getAliasedFile(MODEL));

		storeLoss();
		
		setStatus("Starting applying " + getQualifiedName(receivedConf));
		return applyModel(dtDescriptors.getSliceDescriptorsTable(0, initialSize), receivedConf);
	}

	protected void storeLoss() throws IOException {
	}

	String getPredictionFileName(){
		return PREDICTIONS;
	}

	String[] getCommands(ModelAbstractConfiguration receivedConf, boolean training) throws Exception {
		return new String[] {executablePython, getExeFile(), getAliasedFileName(CFG)};
	}

	@Override
	final protected DataTable applyModelMulti(DescriptorsTable dtDescriptors, MultiLearningAbstractConfiguration receivedConf) throws Exception {

		System.out.println("mols="+dtDescriptors.getDataSize() + " noerr="+dtDescriptors.getRawData().getRowsNoErrorsSize());

		saveModelToFile(receivedConf, getAliasedFileName(MODEL));

		int initialMolSize = dtDescriptors.getDataSize();

		saveConfig(receivedConf, dtDescriptors, false, false);  // order is important (before descriptors) to initialize DeepChem method

		if(saveAggregatedData(getAliasedFileName(APPLYFILE), dtDescriptors, null, receivedConf) == 0)    // nothing can be saved
			return restoreFailedRows(new DataTable(true), dtDescriptors);

		Exception ee = null;

		int maxtrial = noGPU()?1:3;

		if(isRunningTest())maxtrial = 1;

		for(int trial = 0; trial < maxtrial ; trial++) try{
			boolean forceCPU =  trial%2 == 1;

			if(trial > 0)
				setStatus("Applying " + getQualifiedName(receivedConf) + (trial>0?" attempt:" + trial:"") + " using "+(forceCPU?"CPU":"GPU:"+getGPUCard(forceCPU)));

			saveConfig(receivedConf, dtDescriptors, false, forceCPU);

			executablePython = runPythonWithConda(getCommands(receivedConf,false), getPredictionFileName(), getCondaEnvironment());

			DataTable res = readResultValues(getPredictionFileName(), receivedConf, startPositionResults); // only not failed predictions are here

			System.out.println("Augmenting predictions (if any):");
			int columns = res.getColumnsSize();  // number of predictions
			for(int i=0;i<columns;i++)res.addColumn(QSPRConstants.DM + (i == 0 && columns == 1? "": i) + QSPRConstants.STDEVAUG);

			ArrayList<Double> sortarray = new ArrayList<Double>();
			int alreadyAnalysed = dtDescriptors.getRawData().getRowsNoErrorsSize(); // does not contain augmented data, only initial ones

			for(int i = 0, resi = 0; i < initialMolSize; i++)if(!dtDescriptors.getRawRow(i).isError()) // OK we can analyze it
			{
				Integer augmentations  = (Integer)dtDescriptors.getRawRow(i).getAttachment(QSPRConstants.AUGMENTATIONS);
				if(augmentations != null && augmentations > 0) {
					for(int c = 0; c < columns ; c++) {
						sortarray.clear();
						sortarray.add((Double) res.getValue(resi, c));

						for(int j = 0; j < augmentations ; j++) 
							sortarray.add((Double) res.getValue(alreadyAnalysed + j, c));

						double mean = 0, stdev = 0;
						for(Double val: sortarray) mean += val;
						mean /= (augmentations+1);

						for(Double val: sortarray) stdev += (mean - val)*(mean - val);

						res.setValue(resi, c, mean);
						res.setValue(resi, c + columns, Math.sqrt(stdev/augmentations));
					}
					alreadyAnalysed += augmentations;
				}
				resi++;
			}
			if(alreadyAnalysed != res.getRowsSize())throw new IOException("Not correct number of processed records");
			res = res.getSlice(0, dtDescriptors.getRawData().getRowsNoErrorsSize()); // indeed only those were processed

			return restoreFailedRows(res, dtDescriptors);
		}catch(Throwable e) {
			ee= new IOException(e.getMessage());
			e.printStackTrace();
			setStatus("Failed with " + e.getMessage());
		}
		throw ee;
	}

	private DataTable restoreFailedRows(DataTable results, DescriptorsTable dtDescriptors) throws IOException {

		if(results.getRowsSize() != dtDescriptors.getDataSize())
			for(int i=0; i< dtDescriptors.getDataSize();i++)
				if(dtDescriptors.getRawRow(i).isError()) {
					AbstractDataRow row = results.createRow();
					row.setError(dtDescriptors.getRawRow(i).detailedStatus);
					results.insertRow(i, row);
				}
		return results;
	}

	protected String methodType(){
		return supportedTaskType.toUpperCase();
	}

	@Override
	public int saveAggregatedData(String filename, DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration conf) throws Exception {

		NoDescriptors receivedConf = (NoDescriptors) conf;

		if(dtExpValues == null) {
			dtDescriptors = dtDescriptors.getSliceDescriptorsTable(0, dtDescriptors.getDataSize()); // for model training the same size as dtExpValuesO
			prepareSMILES(dtDescriptors, null, receivedConf.isInternalAugmentation()?1:receivedConf.getAugmentApply(), receivedConf);
		}
		else {
			dtDescriptors = dtDescriptors.getSliceDescriptorsTable(0, dtExpValues.getDataSize()); // for model training the same size as dtExpValuesO
			dtExpValues = (LabelsTable) dtExpValues.getSlice(0, dtExpValues.getDataSize());
			prepareSMILES(dtDescriptors, dtExpValues, receivedConf.isInternalAugmentation()?1:Math.abs(receivedConf.getAugementTraining()),  receivedConf);
		}

		if(dtDescriptors.getRawData().getRowsNoErrorsSize() == 0 )return 0; // no records - will return all errors

		if(dtDescriptors.getRawData().getRowsNoErrorsSize() != dtDescriptors.getDataSize()) { // there are errors: filtering them before saving
			DescriptorsTable dt = dtDescriptors.getSliceDescriptorsTable(0, 0);
			LabelsTable de = dtExpValues == null? null: (LabelsTable) dtExpValues.getSlice(0, 0);
			for(int i=0;i<dtDescriptors.getDataSize();i++) 
				if(!dtDescriptors.getRawRow(i).isError())
				{
					dt.getRawData().addRow(dtDescriptors.getRawRow(i));
					if(de != null) de.getRawData().addRow(dtExpValues.getRawRow(i));
				}
			dtDescriptors = dt;
			dtExpValues = de;
		}

		// saving headers and values 
		BufferedWriter bw = getAliasedBufferedWriter(filename);
		saveHeaders(bw, dtDescriptors, dtExpValues == null?null:dtExpValues.getColumnSize(), conf);
		int savedRecords = saveAggregatedDescriptorsAndValues(dtDescriptors, dtExpValues, bw, conf);
		bw.close();

		saveMethodSpecificData(dtDescriptors, dtExpValues, conf);

		return savedRecords;
	}

	void saveMethodSpecificData(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration conf) throws IOException {
	}

	void saveHeaders(BufferedWriter bw, DescriptorsTable dtDescriptors, Integer outputs, ModelAbstractConfiguration conf) throws IOException {
		// saving all values
		bw.write(QSPRConstants.SMILES_FORMAT.toLowerCase());
		for(int i=0;i<dtDescriptors.getDescriptorsSize();i++)
			bw.write("," + QSPRConstants.DESCRIPTOR + i);

		if(outputs != null) {
			if(outputs == 2 && conf.areClassificationData() && conf instanceof SupportsOneOutputOnly) outputs = 1;
			for(int i=0;i<outputs;i++)
				bw.write("," + QSPRConstants.PREDICTION_RESULT_COLUMN + i);
		}
		bw.write("\n");
	}

	Map<String, ArrayList<String>>  augmentData(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, int augmentation, int augments[], NoDescriptors receivedConf) throws Exception {

		boolean stereochemistry = receivedConf instanceof Stereochemistry? ((Stereochemistry)receivedConf).requiresStereochemistry():true;
		boolean training = dtExpValues != null;

		Map<String,ArrayList<String>> smilesArray = new HashMap<String,ArrayList<String>>();

		int initialMol = dtDescriptors.getDataSize();

		int maxprop = 0;
		if (augments!=null) {
			for(int j = 0; j<augments.length;j++) {
				if(maxprop<augments[j])
					maxprop=augments[j]; // selecting property with maximum number of records
				augments[j] *= -1;
			}
		}

		for(int n = 0; n < 1 || (augments != null && n < augments.length); n++) {

			Map<String,Set<String>> smiles = new HashMap<String,Set<String>>();

			int prop = -1;
			if(augments != null) { // select property with maximum number of repetitions
				for(int j = 0; j<augments.length;j++)
					if(augments[j]<0) {
						if(prop == -1) prop = j;
						if(augments[j]>augments[prop])prop=j;
					}
				if(prop == -1) break; // there were properties without any data - we just skip them ...
				augments[prop] = (int)(-maxprop*augmentation/augments[prop]+0.4999); // converting to maximum number of values per property
				if(augments[prop] < 1) augments[prop] = 1;
			}

			BufferedWriter writer = getAliasedBufferedWriter(QSPRConstants.SMILES_FORMAT);

			for(int i=0;i<initialMol;i++) {
				AbstractDataRow row = dtDescriptors.getRawRow(i);
				String smile = (String) row.getAttachment(QSPRConstants.SMILES_ATTACHMENT);

				if(smile ==  null || smile.length() == 0 ||  row.isError() || smile.toUpperCase().contains(QSPRConstants.ERROR) || smilesArray.containsKey(smile) || 
						(augments != null && dtExpValues.isMissedValue(i,prop)))continue;

				smiles.put(smile, new LinkedHashSet<String>());
				smilesArray.put(smile, new ArrayList<String>());
				writer.append(smile+"\n");
			}

			writer.flush();
			writer.close();

			if(smiles.size() == 0) continue;
			System.out.println("using augmentation for " + smiles.size());

			int augm = augmentation;
			if(augments != null) {
				augm = augments[prop];
				System.out.println("using augmentation " + augm + " for property " + prop);
			}

			String[] commands = new String[] {null, getFilter(), "--infile", QSPRConstants.SMILES_FORMAT, "--outfile",output,"--augment",""+(augm>1?2*augm:1),"--isomeric",
					stereochemistry?"True":"False","--methodtype",methodType()};

			if(!receivedConf.isSanitize()) { // by default sanitization is done; thus we need to disable it
				commands = OCHEMUtils.appendString(commands, "-n");
				commands = OCHEMUtils.appendString(commands, "True");
			}

			runPythonWithConda(commands, output, CONDA.RDKIT);

			BufferedReader br = getAliasedBufferedReader(output);

			int goodSmiles = 0;

			String line;
			while ((line = br.readLine()) != null) {
				if(line.toUpperCase().contains(QSPRConstants.ERROR))continue;
				String pieces[] = line.split(",");
				if(pieces.length != 2)continue;
				String augSmile = pieces[1];
				if(receivedConf.isPeptide(pieces[0]) != null)augSmile = receivedConf.isPeptide(pieces[0]);
				if(receivedConf.isGoodSMILES(augSmile, training) != null)continue;
				if(!smiles.containsKey(pieces[0]))continue;
				Set<String> s = smiles.get(pieces[0]);
				s.add(augSmile);
				goodSmiles++;
			}
			br.close();


			if(training && goodSmiles == 0)
				throw new CriticalException("Calculation of all molecules failed");

			if(! ((ModelAbstractConfiguration)receivedConf).skipSizeCheck() &&
					(training && goodSmiles<((ModelAbstractConfiguration)receivedConf).requireMinimumRecords()/2 && !isRunningTest()))
				throw new CriticalException(QSPRConstants.IMMEDIATE_FAILURE + "Most of molecules failed to be processed and available number \"" +  goodSmiles + "\" is less than required by this method: " + ((ModelAbstractConfiguration)receivedConf).requireMinimumRecords());

			for(String aug:smiles.keySet()){
				if(smiles.get(aug).size()==0)continue;
				ArrayList<String> list = smilesArray.get(aug);
				for(String smile:smiles.get(aug)){
					list.add(smile);
					if(list.size()==augm)break;
				}
				if(augments != null) // we need to do balancing -- add duplicates
					for(int inisize = list.size();list.size()!=augm;)
						list.add(list.get(ThreadLocalRandom.current().nextInt(0, inisize)));
			}
		}

		return smilesArray;

	}

	void prepareSMILES(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, int augmentation, NoDescriptors receivedConf) throws Exception {

		if(dtExpValues != null && dtExpValues.getDataSize() != dtDescriptors.getDataSize())
			throw new IOException("sites of data with and without labels should be identicall!");

		int initialMol = dtDescriptors.getDataSize();
		String originalMols[] = new String[initialMol]; // strange, but datatable somehow contains reference to attachment, which is updated in several rows simultaneously.. Really strange...
		int mol = 0;
		for(int i=0;i<initialMol;i++) { 
			originalMols[i] = (String) dtDescriptors.getRawRow(i).getAttachment(QSPRConstants.SMILES_ATTACHMENT);
			if(originalMols[i] != null)mol++;
		}

		if(dtExpValues != null && mol < ((ModelAbstractConfiguration)receivedConf).requireMinimumRecords()/2 && !(isRunningTest()))
			throw new CriticalException("Number of molecules with non empty SMILES " + mol + " < " + ((ModelAbstractConfiguration)receivedConf).requireMinimumRecords());

		int augments[] = receivedConf.getBalanceData()&& dtExpValues != null?dtExpValues.countNonMissedNotEmptyClassValues():null;

		Map<String,ArrayList<String>>  smiles = augmentData(dtDescriptors, dtExpValues, augmentation, augments, receivedConf);

		double all = 0;
		if(augments != null) {
			for(int i=0;i<initialMol;i++) { // counting the number of data to be saved for each property
				AbstractDataRow row = dtDescriptors.getRawRow(i);
				String smile = originalMols[i];
				ArrayList<String> augsmiles = smiles.get(smile);
				if(augsmiles == null || augsmiles.isEmpty()) { // all conversion failed -- skipping this row
					String mess = receivedConf.isGoodSMILES(smile,dtExpValues != null);
					row.setError(mess == null?"RDKit failed to process: " + smile +" You can try to avoid this error by enabling: \"Skip sanitisation\" option":mess);
					continue;
				}

				int property = dtExpValues.whichProperty(i);
				int num = augsmiles.size();
				num = num > augments[property]?augments[property]:num;
				all += num;
			}

			double coef = initialMol*augmentation / all;

			for(int i=0;i<augments.length;i++) {
				int previous = augments[i];
				augments[i] = (int)(augments[i]*coef+0.4999);
				augments[i] = augments[i] < 1? 1:augments[i];
				System.out.println("ajusted for " + i + " " + previous + " --> " + augments[i]);
			}
		}

		for(int i=0;i<initialMol;i++) {
			AbstractDataRow row = dtDescriptors.getRawRow(i);
			String smile = originalMols[i];
			ArrayList<String> augsmiles = smiles.get(smile);
			if(augsmiles == null || augsmiles.isEmpty()) { // all conversion failed -- skipping this row
				row.setError("RDKit failed to process: " + smile);
				continue;
			}

			row.addAttachment(QSPRConstants.SMILES_ATTACHMENT, augsmiles.get(0)); // original row, not changing it - just update SMILES

			int num = augsmiles.size();
			if(augments != null) { // for over-represented properties we limit number of records
				int property = dtExpValues.whichProperty(i);
				num = num > augments[property]?augments[property]:num;
				Collections.shuffle(augsmiles);
			}

			for(int n = 1; n < num; n++) {
				AbstractDataRow rowNew = row.getDeeperCopy();
				rowNew.addAttachment(QSPRConstants.SMILES_ATTACHMENT, augsmiles.get(n));
				dtDescriptors.getRawData().addRow(rowNew);
				if(dtExpValues != null)
					dtExpValues.getRawData().addRow(dtExpValues.getRawData().getRow(i));
			}
			row.addAttachment(QSPRConstants.AUGMENTATIONS, augsmiles.size()-1); // number of augmentation
		}
		all = dtDescriptors.getRawData().getRowsSize();
		int good = dtDescriptors.getRawData().getRowsNoErrorsSize();
		System.out.println("Augmented SMILES: total " + all + " good: " + good + " failed: " + 
				(dtDescriptors.getDataSize() - good) + " Average = " + (all != 0?Math.round(((float)good)/initialMol +0.4999):"0")+ " for initial: " + initialMol);

	}
}
