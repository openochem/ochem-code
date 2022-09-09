
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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

import com.eadmet.exceptions.CriticalException;

import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.KGCNNConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.utils.QSPRConstants;

/**
 * Requires specific handling, since instead of descriptors SDF data are
 * provided
 * 
 * @author itetko
 *
 */

public class KGCNNServer extends SmilesOnlyAbstractServer {

	KGCNNConfiguration config;
	final static String LOSS = "loss.csv";

	/** Initializes supported task and workflow. */
	public KGCNNServer() {
		supportedTaskType = QSPRConstants.KGCNN;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
		startPositionResults = 0;
		batchApplySize = 10000;
	}

	@Override
	public boolean isCritical(String message) {
		if(message.contains("NotImplementedError"))
			throw new CriticalException("NotImplementedError: " + message);
		return super.isCritical(message);
	}


	@Override
	protected String[] getMessages() {
		return new String[] { "EPOCH:"," loss:","val_loss:"};
	}

	@Override
	protected CONDA getCondaEnvironment(){
		return CONDA.RDKIT;
	}

	@Override
	protected String getGPUCard(boolean forceCPU) { // always card 0, since environment is used to restrict only one card, which becomes 0
		return "" + (noGPU() || forceCPU ? NO_GPU : 0);
	}


	@Override
	protected void storeLoss() throws IOException {
		File f = new File(getAliasedFileName(LOSS));
		if(f.exists()) {
			BufferedReader br = exeRunner.getAliasedBufferedReader(getAliasedFileName(LOSS));

			String mes[] = getMessages();

			String line = br.readLine();
			while ((line = br.readLine()) != null)
			{
				line = line.trim();
				if(line.length() < 10)continue;
				String[] predictions = line.split(",");
				scanStatus(mes[0]+" " + predictions[0]+
						" "+mes[1]+" " + predictions[1]+
						" "+ mes[2]+" " + predictions[3]);
			}
			br.close();

		}
	}

	@Override
	void saveHeaders(BufferedWriter bw, DescriptorsTable dtDescriptors, Integer outputs, ModelAbstractConfiguration conf) throws IOException {
		KGCNNConfiguration config = (KGCNNConfiguration) conf;
		super.saveHeaders(bw, dtDescriptors, outputs(config), conf);
	}

	int outputs(KGCNNConfiguration conf) {
		boolean isRegression = conf.containAlsoRegressionData();

		int outputs = conf.OutputValues();
		if(!isRegression && outputs == 2)
			return 1;

		return outputs;
	}

	@Override
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean train, boolean forceCPU) throws IOException{

		KGCNNConfiguration conf = (KGCNNConfiguration) configuration;

		config = conf;
		
		config.sanitize = true;

		File f = new File(getAliasedFileName(CFG));
		if(f.exists()) {
			File ff = new File(getAliasedFileName(CFG+".train"));
			if(!ff.exists())
				f.renameTo(ff);
		}

		boolean isRegression = conf.containAlsoRegressionData();

		String configur=
				"[Task]\n"
						+ "train_mode = " + (train? "True":"False") +"\n"
						+ "model_file = model.tar\n"
						+ (train?"train_data_file = " + DATAFILE : "apply_data_file = "+ APPLYFILE) + "\n"
						+ "result_file= " + PREDICTIONS + "\n"
						+ "\n"
						+ "[Details]\n"
						+ "batch = " + conf.batch +"\n"
						+ "architecture_name = " + conf.method +"\n"
						+ "nbepochs = " + conf.nepochs+ "\n"
						+ "gpu = " + getGPUCard(forceCPU) + "\n"
						+ "seed = " + conf.getSeed() + "\n"
						+ "output_dim = " + conf.OutputValues() + "\n" 
						//+ "activation = sigmoid\n"
						//+ "classloss = BCEmask\n" 
						+ "overwrite = True\n"
						+ "activation = " + (isRegression?"linear":"sigmoid") + " \n"
						+ "classloss=" + (isRegression?"RMSEmask":"BCEmask")+ "\n"
						+ "classification=" + (isRegression?"False":"True") + "\n"
						;


		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		writer.write(configur );
		writer.close();
	}


}
