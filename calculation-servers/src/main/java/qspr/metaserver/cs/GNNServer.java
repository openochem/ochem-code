
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
import java.io.File;
import java.io.IOException;

import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.GNNConfiguration;
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

public class GNNServer extends SmilesOnlyAbstractServer {
	/** Initializes supported task and workflow. */
	public GNNServer() {
		supportedTaskType = QSPRConstants.GNN;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
		MODEL = "model";
	}

	@Override
	protected String[] getMessages() {
		return new String[] { "Epoch","normalized train loss","val loss"};
	}

	@Override
	protected String getGPUCard(boolean forceCPU) { // always card 0, since environment is used to restrict only one card, which becomes 0
		return "" + (noGPU() || forceCPU ? NO_GPU : 0);
	}

	@Override
	protected CONDA getCondaEnvironment(){
		return CONDA.RDKIT;
	}

	@Override
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean train, boolean forceCPU) throws IOException{

		GNNConfiguration conf = (GNNConfiguration) configuration;

		File f = new File(getAliasedFileName(CFG));
		if(f.exists()) {
			File ff = new File(getAliasedFileName(CFG+".train"));
			if(!ff.exists())
				f.renameTo(ff);
		}

		BufferedWriter writer = getAliasedBufferedWriter(CFG);

		writer.write("[parameters]" +
				"\ntrain = " + (train? "True":"False") +
				"\ndim = " + conf.dim + 
				"\nbatch = " + conf.batch +
				"\npatience =  "  + conf.patience+
				"\npatience_early = " + conf.patienceearly +
				"\nlr_decay = " + conf.lr_decay +
				"\nn_features = 48" +
				"\nn_epochs = " + conf.nepochs+
				"\nnum_atoms = True" +
				"\nuse_hydrogen_bonding = False" +
				"\nneighbour_dist = False" +
				"\nuse_acid_base = False" +
				"\npartial_charge = False" +
				"\nn_iterations = 10" +
				"\nmodel = " + conf.gnn
				);

		int outputs = conf.OutputValues();

		writer.write("\ntask_names = ");
		for(int i = 0; i< outputs; i++)
			writer.write(QSPRConstants.PREDICTION_RESULT_COLUMN + i + ((i != (outputs-1))?",":""));

		writer.write(
				"\ngpu = "+ getGPUCard(forceCPU)  +
				"\nseed = " +conf.getSeed() +"\n"
				);

		writer.close();
		return;
	}

}
