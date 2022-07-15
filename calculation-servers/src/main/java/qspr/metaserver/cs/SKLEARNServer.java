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

import java.io.BufferedWriter;
import java.io.IOException;

import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.SKLConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;


/**
 *  Is in incorrect place! This is not multilearning server
 * @author itetko
 *
 */

public class SKLEARNServer extends MachineLearningExecutableAbstractServer
{
	final String MODEL =  "best_model.hdf"; 
	final String DATAFILE = "train.libsvm"; 
	final String APPLYFILE = "apply.libsvm"; 
	final String PREDICTIONS = "apply_results.txt";

	final String applyModel = "apply.py";
	final String trainModel = "train.py";
	final String CFG = "model.cfg";

	/** Initializes supported task and workflow. */
	public SKLEARNServer()
	{
		supportedTaskType = QSPRConstants.SKLEARN;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}

	String[] prepareCommands(ModelAbstractConfiguration receivedConf, boolean train) throws IOException{
		saveConfig((SKLConfiguration)receivedConf,train);
		String[] commands = { null, "run.py", CFG};
		return commands;
	}

	void saveConfig(SKLConfiguration conf, boolean traningMode) throws IOException{
		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		writer.write("[Task]\n");
		writer.write("\nbase_model =  " + ("" + conf.method).toLowerCase());
		writer.write("\ntrain_mode = " + traningMode);
		writer.write("\nmodel_file = " + MODEL);
		writer.write("\ntask_names = ");
		for(int i = 0; i< conf.OutputValues() || i == 0; i++)
			writer.write(QSPRConstants.PREDICTION_RESULT_COLUMN + i + ((i != (conf.OutputValues() -1))?",":""));
		writer.write("\ntrain_data_file = " + DATAFILE);
		writer.write("\napply_data_file = " + APPLYFILE);
		writer.write("\nresult_file = " + PREDICTIONS);

		writer.write("\n\n[Details]");
		writer.write("\nseed = " + conf.getSeed());
		writer.write("\nparallelize = False");

		writer.write("\n");
		writer.close();
	}

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues,
			ModelAbstractConfiguration receivedConf) throws Exception {
		
		saveDataSVM(dtDescriptors, dtExpValues, DATAFILE);
		saveDataSVM(dtDescriptors, null, APPLYFILE);

		runPythonWithConda(prepareCommands(receivedConf, true), MODEL,null);
		receivedConf.storeModel(getAliasedFile(MODEL));
		return applyModel(dtDescriptors, receivedConf);
	}

	@Override
	protected DataTable applyModelBatch(DescriptorsTable dtDescriptors, ModelAbstractConfiguration receivedConf) throws Exception {
		saveDataSVM(dtDescriptors, null, APPLYFILE);
		saveModelToFile(receivedConf, getAliasedFileName(MODEL));
		runPythonWithConda(prepareCommands(receivedConf, false), PREDICTIONS, null);
		return readResultValuesSingle(PREDICTIONS,receivedConf.optionsNumber,0);
	}

	private void saveDataSVM(DescriptorsTable original, LabelsTable vals, String file) throws IOException, InterruptedException{
		BufferedWriter bw  = getAliasedBufferedWriter(file);
		int size = vals == null? original.getDataSize() : vals.getDataSize();
		for(int i =0; i < size ; i++)
			original.saveOneRowSVM(i, vals == null ? null: new String[]{vals.getOneValueOnlyString(i)}, bw);
		bw.close();
	}
}
