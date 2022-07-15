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
import java.io.File;
import java.io.IOException;

import com.eadmet.exceptions.CriticalException;

import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.ATTFPConfiguration;
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

public class ATTFPServer extends SmilesOnlyAbstractServer {
	/** Initializes supported task and workflow. */
	public ATTFPServer() {
		supportedTaskType = QSPRConstants.ATTFP;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
		startPositionResults = 1;
		batchApplySize = 10000;
	}

	@Override
	protected CONDA getCondaEnvironment(){
		return CONDA.RDKIT;
	}

	@Override
	public boolean isCritical(String message) {
		if(message.contains("cannot reindex on an axis with duplicate labels frame"))
			throw new CriticalException("cannot reindex on an axis with duplicate labels frame: " + message);
		return super.isCritical(message);
	}

	@Override
	protected String getGPUCard(boolean forceCPU) { // always card 0, since environment is used to restrict only one card, which becomes 0
		return "" + (noGPU() || forceCPU ? NO_GPU : 0);
	}

	@Override
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean train, boolean forceCPU) throws IOException{

		ATTFPConfiguration conf = (ATTFPConfiguration) configuration;

		conf.sanitize = true;
		
		File f = new File(getAliasedFileName(CFG));
		if(f.exists()) {
			File ff = new File(getAliasedFileName(CFG+".train"));
			if(!ff.exists())
				f.renameTo(ff);
		}

		int outputs = configuration.OutputValues();
		if(outputs == 2 && conf.areClassificationData()) outputs = 1;

		int batch = conf.batch;

		if(!train && batch > 32) batch = 32; 

		String configur=
				"[Task]\n"
						+ "train_mode = " + (train? "True":"False") +"\n"
						+ "model_file = model.tar\n"
						+ (train?"train_data_file = " + DATAFILE : "apply_data_file = "+ APPLYFILE) + "\n"
						+ "result_file= " + PREDICTIONS + "\n"
						+ "\n"
						+ "[Details]\n"
						+ "model = AttentionFP-v2\n"
						+ "batch = " + batch +"\n"
						+ "lr = "+ conf.lr +"\n"
						+ "nbepochs = " + conf.nepochs+ "\n"
						+ "patience_reduce = " + conf.patience_reduce+ "\n"
						+ "patience_early = " + conf.patience_early+ "\n"
						+ "dropout = " + conf.dropout+ "\n"
						+ "fp_dim = " + conf.fp_dim+  "\n"
						+ "cosineT = " + conf.cosineT+  "\n"
						+ "radius= " + conf.radius +" \n"
						+ "T= " + conf.T +" \n"
						+ "weight_decay= " + conf.weight_decay +" \n"

						+ "output_dim= " + outputs +" \n"

						+ "cosine = "  + (conf.cosine?"True":"False")+ "\n"
						+ "early = "  + (conf.early?"True":"False")+ "\n"

						+ "gpu = " + getGPUCard(forceCPU) + "\n"
						+ "simpleO =" + (conf.simpleO !=null && conf.simpleO?"True":"False") + "\n"
						+ "singleT =" + (conf.singleT !=null && conf.singleT?"True":"False") + "\n"
						+ "lngru =" + (conf.lngru !=null && conf.lngru?"True":"False") + "\n"
						+ "isRAdam = True\n"
						+ "best = True\n" 
						;


		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		writer.write(configur );
		writer.close();
		return;
	}

}
