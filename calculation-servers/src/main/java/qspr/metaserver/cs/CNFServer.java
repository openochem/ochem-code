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

import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.configurations.CNFConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.workflow.utils.QSPRConstants;

public class CNFServer extends SmilesOnlyAbstractServer{
	@Override
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean train, boolean forceCPU) throws IOException{

		CNFConfiguration conf = (CNFConfiguration) configuration;

		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		if(train) {
			writer.write("[Task]"+
					"\ntrain_mode = True" +
					"\nmodel_file = "+ MODEL +
					"\ntrain_data_file = " + DATAFILE
					);

			int outputs = conf.OutputValues();


			//if(conf.getImplicit() != null)
			//	outputs = conf.OutputValues()/2;

			writer.write("\ntask_names = ");
			for(int i = 0; i< outputs; i++)
				writer.write(QSPRConstants.PREDICTION_RESULT_COLUMN + i + ((i != (outputs-1))?",":""));

			if(conf.labelWeighting != null)
				writer.write("\nweights = " + WEIGHTS);

			saveCommonCfg(writer, conf, getGPUCard(forceCPU));
		}
		else
			writer.write("[Task]\n"+
					"train_mode = False\n"+
					"apply_data_file = " + APPLYFILE +"\n" +
					"result_file = " + PREDICTIONS + "\n" +
					"model_file = "+ MODEL + "\n" +
					"\n[Details]\n"+
					"gpu = "+ getGPUCard(forceCPU) + "\n"
					);

		writer.write("\n");
		writer.close();
		return;
	}


	/** Initializes supported task and workflow. */
	public CNFServer()
	{
		supportedTaskType = QSPRConstants.CNF;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
		batchApplySize = 10000;
	}


	public static void saveCommonCfg(BufferedWriter writer, CNFConfiguration conf, String gpu) throws IOException {
		String error = conf.containAlsoRegressionData()?"rmse":"entropy";
		writer.write("\n\n[Details]\n"
				+ "gpu = "+ gpu  + "\n" 
				+ "seed = "+ conf.getSeed() + "\n"
				+ "n_epochs = " + conf.nepochs + "\n"
				+ "fp_layer_dim = " + conf.FPDim + "\n" // FPDim
				+ "n_fp_layers = " +conf.nLayers +"\n" //nLayers
				+ "mlp_dims = " + conf.FFNET_dim+"\n" // FFNET_dim
				+ "batch_size = " + conf.batch +"\n"
				+ "learning_rate = " + conf.rate +"\n" 
				+ "early = " + conf.early +"\n"
				+ "error_function = " + error + "\n"
				+ "augment = " + (conf.augmentation == null || conf.augmentation != -1 ?0: -1)  + " \n"
				+ "tokenizer_type  = " + conf.getTokeniser().toString().toLowerCase() + "\n"
				+ "dropout  = " + conf.getDropout()+ "\n"
				+ "relu_type  = " + conf.getActivationFunction()+ "\n"
				//+ "canonize = False\n"
				// + (conf.highway != null ? "highway = " + conf.highway + "\n" : "")
				// + ("allow_overfit  = " + (conf.shuffleAllData()?"False":"True") + "\n")
				);
	}


}
