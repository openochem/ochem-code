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
import java.util.LinkedHashMap;

import com.eadmet.utils.FileUtils;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.OSType;

import qspr.metaserver.configurations.KPLSConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;


public class KPLSServer extends IgorsAbstractServer
{	
	protected int getResColumn() {
		return 2;
	}

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws Exception
	{
		KPLSConfiguration cfg = (KPLSConfiguration) receivedConf;
		cfg.latentVariables  = dtDescriptors.getColumnSize()>cfg.latentVariables?cfg.latentVariables:dtDescriptors.getColumnSize();

		FileUtils.saveStringToFile("  "+cfg.latentVariables,getAliasedFileName("num_eg.txt"));

		//create data file for KNN Training
		BufferedWriter writer = getAliasedBufferedWriter(DATAFILE);
		writeSimpleIgorFormatData(dtDescriptors, dtExpValues, writer);
		writer.close();

		String[] commands = new String[] { getExeFile(), DATAFILE, "--SCALE_PAT"};
		executeBinary(commands, DATAFILE+".txt", 0);
		if(cfg.latentVariables == 0) {
			setStatus("Starting LAMBDA_TUNE ...");
			commands[1] = DATAFILE+".txt";  commands[2] = cfg.useLambda != null && cfg.useLambda? "--LAMBDA_TUNE2" :"--LAMBDA_TUNE";
			executeBinary(commands, "num_eg.txt", 0);
			String num = FileUtils.getFileAsString(getAliasedFileName("num_eg.txt")).trim();
			cfg.latentVariables = Integer.valueOf(num);
		}
		setStatus("Starting SIGMA_TUNE ...");
		commands[1] = DATAFILE+".txt"; commands[2] = "--SIGMA_TUNE";
		executeBinary(commands, "sigmas.txt", 0);
		commands[2] = DATAFILE+".txt"; commands[2] = "--TRAIN_KPLS";
		executeBinary(commands, "resultss.xxx", 0);
		commands[1] = "resultss.xxx"; commands[2] = "--DESCALE";
		executeBinary(commands, "results_H.ttt", 0);

		commands[0] = "tar";
		commands[1] = "-cf";
		commands[2] = MODEL;
		commands = OCHEMUtils.append(commands, "ccmatrixx.txt");
		commands = OCHEMUtils.append(commands, "bbmatrixx.txt");
		commands = OCHEMUtils.append(commands, "sigma.dbd");
		commands = OCHEMUtils.append(commands, "resultss.xxx");

		executeBinary(commands, MODEL, 0);

		receivedConf.storeModel(getAliasedFile(MODEL));

		setStatus("KPLS has finished applying the model");
		return applyModel(dtDescriptors, receivedConf);
	}


	@Override
	protected DataTable applyModelBatch(DescriptorsTable dtDescriptors,  ModelAbstractConfiguration receivedConf) throws Exception
	{

		KPLSConfiguration cfg = (KPLSConfiguration) receivedConf;

		BufferedWriter writer = getAliasedBufferedWriter(APPLYFILE);
		writeSimpleIgorFormatData(dtDescriptors, null, 1, writer,null);
		writer.close();

		if(cfg.getTrainingSetDescriptors() != null) { // if we apply model immediately after the training
			DescriptorsTable dtTraining = new DescriptorsTable(cfg.getTrainingSetDescriptors(),cfg,0);
			LabelsTable dtExpValues = new LabelsTable(cfg.getTrainingSetValues(), cfg);

			writer = getAliasedBufferedWriter(DATAFILE);
			writeSimpleIgorFormatData(dtTraining, dtExpValues, writer);
			writer.close();

			String[] datas = new String[] { getExeFile(), DATAFILE, "--SCALE_PAT"};
			executeBinary(datas, DATAFILE+".txt", 0);
		}

		saveModelToFile(receivedConf, getAliasedFileName(MODEL));
		String untar[]=new String[] { "tar","-xf",MODEL};
		executeBinary(untar, "ccmatrixx.txt", 0);

		String[] commands = new String[] { getExeFile(), APPLYFILE, "--SCALE_TES"};
		executeBinary(commands, APPLYFILE+".txt", 0);
		commands[1] = APPLYFILE+".txt"; commands[2] = "--TEST_KPLS";
		executeBinary(commands, "resultss.ttt", 0);
		commands[1] = "resultss.ttt"; commands[2] = "--DESCALE";
		executeBinary(commands, "results.ttt", 0);

		return readResult("results.ttt");
	}


	@Override
	protected void writeSimpleIgorFormatData(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, int outputs, BufferedWriter writer, LinkedHashMap<String,Integer> columnsToKeep)
			throws IOException {
		int allMolecules = dtExpValues == null? dtDescriptors.getDataSize(): dtExpValues.getDataSize();

		VALUESSEPARATOR =",  ";

		boolean cols[]=new boolean[dtDescriptors.getDescriptorsSize()];

		for(int i=0;i<cols.length;i++)
			cols[i]=columnsToKeep==null?true:columnsToKeep.containsKey(dtDescriptors.columnNames.get(i));

		// saving molecules one by one
		for (int mol = 0; mol < allMolecules; mol++){
			writer.append("  "+mol+VALUESSEPARATOR);

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
			dtResult.addRow();
			String[] prediction = line.split("\\s+");
			dtResult.setValue(QSPRConstants.PREDICTION_RESULT_COLUMN, Double.valueOf(prediction[getResColumn()]));
		}
		graphics.close();

		return dtResult;

	}


	public KPLSServer(){
		supportedTaskType = QSPRConstants.KPLS;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}
}
