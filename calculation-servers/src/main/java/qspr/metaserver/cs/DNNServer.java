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

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.net.InetAddress;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.utils.OCHEMUtils;

import qspr.metaserver.configurations.DNNConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

public class DNNServer extends MultiLearningAbstractServer{

	private static final int TRIALS = 5; // to prevent crashing of DNN on linux servers
	private static final int LARGETRIALS = 50;
	final String MODEL =  "best_model.hdf"; 
	final String DATAFILE = "train.libsvm"; 
	final String DATALOG = "train.log"; 
	final String APPLYFILE = "apply.libsvm"; 
	final String PREDICTIONS = "apply_results.txt";

	final String applyModel = "apply.py";
	final String trainModel = "train.py";

	String[] prepareCommands(DNNConfiguration configuration, boolean train) throws IOException{

		String[] commands; String python =null;

		if(train) {

			commands = new String[] {python, getAliasedFileName(trainModel), "--batch_size",
					""+configuration.batchSize, "--clear_cupy_cache", "--internal_train_test_ratio",""+configuration.ratio,
					"--optimizer",configuration.optimizertype,"--seed",""+ configuration.getSeed(), "--gpu_card", getGPUCard(false),
					"--n_epochs",""+configuration.epochs
			};

			if(configuration.areClassificationData()) {
				commands = OCHEMUtils.append(commands, "--lossfunc");
				commands = OCHEMUtils.append(commands, "sce");
			}

			commands = addOptional(commands, configuration.additionalParam );

		}else 
			commands = new String[] {python, getAliasedFileName(applyModel), "--gpu_card", getGPUCard(false)};

		if(noGPU() || gpuCard == NO_GPU)
			commands = OCHEMUtils.append(commands,"--use_cpu");

		commands = OCHEMUtils.append(commands, configuration.modeltype);

		if(configuration.modeltype.equals("dense_exp")){
			commands = OCHEMUtils.append(commands, "--depth");
			commands = OCHEMUtils.append(commands, "5");
			commands = OCHEMUtils.append(commands, "--alpha");
			commands = OCHEMUtils.append(commands, "2");
		}

		commands = addOptional(commands, configuration.additionalParamModel);

		return commands;
	}

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors,
			LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf)
					throws Exception {
		DescriptorsTable dtOriginal = dtDescriptors;

		DNNConfiguration configuration = (DNNConfiguration) receivedConf;

		saveAggregatedData(getAliasedFileName(DATAFILE), dtDescriptors, dtExpValues, receivedConf);
		saveAggregatedData(getAliasedFileName(APPLYFILE),dtDescriptors, null, receivedConf);

		String [] commands = prepareCommands(configuration, true);

		int trials = TRIALS;
		for(int trial=0; trial <= trials; trial++) try{
			if(trial > 0) setStatus("Trying to calculate model - attempt " + trial + " out of " + trials);
			runPythonWithConda(commands, MODEL,CONDA.RDKIT);
			File f = getAliasedFile(DATALOG);
			if(!f.exists())throw new IOException("log file is absent; program has crashed");
			break;
		}catch(CriticalException e){
			badServer(e, trial,trials);
		}
		catch(Exception e){
			badServer(e, trial,trials);
		}

		receivedConf.storeModel(getAliasedFile(MODEL));

		return applyModel(dtOriginal, receivedConf);
	}

	void badServer(Exception e, int trial, int trials) throws Exception {
		setStatus("There was an exception for trial " + trial + " out of " + trials);
		if(trial > TRIALS)Thread.sleep(1000*trial);
		if(trial == trials)throw e;
	}

	@Override
	protected DataTable applyModelMulti(DescriptorsTable dtDescriptors, MultiLearningAbstractConfiguration receivedConf) throws Exception {
		DataTable res = null;

		saveAggregatedData(getAliasedFileName(APPLYFILE),dtDescriptors,null,receivedConf);
		saveModelToFile(receivedConf, getAliasedFileName(MODEL));

		String [] commands = prepareCommands((DNNConfiguration) receivedConf, false);

		int trials = InetAddress.getLocalHost().getHostName().contains("dc01cmsc2241")?LARGETRIALS:TRIALS;

		for(int trial = 0; trial <= trials; trial++) try{
			if(trial > 0) setStatus("Trying to apply model - attempt " + trial + " out of " + trials);
			runPythonWithConda(commands, PREDICTIONS,CONDA.RDKIT);
			res = readResultValues(PREDICTIONS,receivedConf);
			break;
		}catch(CriticalException e){
			badServer(e, trial,trials);
		}
		catch(Exception e){
			badServer(e, trial,trials);
		}

		return res;
	}

	@Override
	protected void saveOneRow(int mol, DescriptorsTable dtDescriptors, String[] values, Writer bw) throws IOException {
		dtDescriptors.saveOneRowSVM(mol, values, bw);
	}

	private String[] addOptional(String[] commands,
			String additionalParam) {

		if(additionalParam != null)
			for(String par : additionalParam.split("\\s+"))
				commands = OCHEMUtils.append(commands,par);

		return commands;
	}


	public DNNServer()
	{
		supportedTaskType = QSPRConstants.DNN;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
		batchApplySize = QSPRConstants.MODEL_REPOST_SIZE;
	}

}