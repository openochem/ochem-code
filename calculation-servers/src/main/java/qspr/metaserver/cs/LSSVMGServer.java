
/**
 * Neural Network Fingerprints
 */

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
import java.io.Writer;

import com.eadmet.utils.OSType;

import qspr.metaserver.configurations.Iterations;
import qspr.metaserver.configurations.LSSVMGConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;


/**
 *  Requires specific handling, since instead of descriptors SDF data are provided
 * @author itetko
 *
 */

public class LSSVMGServer extends MultiLearningAbstractServer
{
	static String CFG = "model.cfg";
	static String DATAFILE = "train.csv";
	static String APPLYFILE = "apply.csv";
	static String MODEL = "model.pkl";
	static String PREDICTIONS = "results.csv";
	private static String executablePython = null;
	
	int iterations = 0;

	String savedStatus = "";

	/** Initializes supported task and workflow. */
	public LSSVMGServer()
	{
		supportedTaskType = QSPRConstants.LSSVMG;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);	
		batchApplySize = 5000;
	}

	@Override
	protected String[] getMessages() {
		return new String[] { "epoch:","train score:","train score:"};
	}

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors,
			LabelsTable dtExpValues, ModelAbstractConfiguration lssvmgConf)
					throws Exception {

		setStatus("Starting calculate " + supportedTaskType);

		saveAggregatedData(getAliasedFileName(DATAFILE), dtDescriptors, dtExpValues, lssvmgConf);
		saveAggregatedData(getAliasedFileName(APPLYFILE), dtDescriptors, null, lssvmgConf);

		saveConfig(lssvmgConf, true, dtExpValues.getColumnSize());
		String[] commands = new String[] {null, getExeFile(), CFG}; // otherwise link to python is not created
		executablePython = runPythonWithConda(commands, MODEL,CONDA.RDKIT);
		lssvmgConf.storeModel(getAliasedFile(MODEL));
		((Iterations)lssvmgConf).setIterations(iterations);

		return applyModel(dtDescriptors,lssvmgConf);
	}

	void saveConfig(ModelAbstractConfiguration configuration, boolean train, int outputs) throws IOException{

		LSSVMGConfiguration conf = (LSSVMGConfiguration) configuration;

		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		if(train) {

			writer.write("[Task]\n"+
					"train_mode = True\n" +
					"model_file = "+ MODEL + "\n" +
					"train_data_file = " + DATAFILE + "\n" +
					"ntargets = " + outputs + "\n" +
					"\n[Details]\n"+
					"gpu = "+ getGPUCard(false) + "\n"+
					"seed = " + conf.getSeed() + "\n"+
					"center = true\n" +
					"pca = false\n" +
					"nfold = " + conf.cv + "\n"+
					"glob_opt = " + (conf.isGlobal()?"True" : "False") + "\n"+
					"metric = rmse\n"+
					"kernel = " + conf.kernel+"\n\n"
					);

			if(conf.additionalParam != null) {
				for (String param : conf.additionalParam.split(","))
					writer.write(param + "\n");
				writer.write("\n");
			}
		}
		else
			writer.write("[Task]\n\n"+
					"model_file = "+ MODEL + "\n" +
					"apply_data_file = " + APPLYFILE +"\n" +
					"gpu = "+ getGPUCard(false) + "\n"+
					"result_file = " + PREDICTIONS + "\n\n" +
					"train_mode = False\n"
					);

		writer.close();
		return;
	}

	@Override
	protected DataTable applyModelMulti(DescriptorsTable dtDescriptors, MultiLearningAbstractConfiguration receivedConf)
			throws Exception {

		saveAggregatedData(getAliasedFileName(APPLYFILE), dtDescriptors, null, receivedConf);
		saveConfig(receivedConf, false, 0);
		saveModelToFile(receivedConf, getAliasedFileName(MODEL));
		String[] commands = new String[] {executablePython, getExeFile(), CFG}; 
		executablePython = runPythonWithConda(commands, PREDICTIONS,CONDA.RDKIT);
		return readResultValues(PREDICTIONS, receivedConf);
	}

	@Override
	protected void saveOneRow(int mol, DescriptorsTable dtDescriptors, String[] values, Writer writer)
			throws IOException {
		String []desc = dtDescriptors.getScaledValuesString(mol);
		for(int i=0;i<desc.length;i++) 
			writer.append(desc[i]+(i!=desc.length-1?",":""));

		if(values != null)
			for(int i=0; i <values.length; i++)
				writer.append("," + values[i]);
		writer.append(OSType.endLine());	
	}

}
