
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

import java.io.IOException;

import com.eadmet.utils.OCHEMUtils;

import qspr.metaserver.configurations.ChemPropConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.workflow.utils.QSPRConstants;


/**
 *  Requires specific handling, since instead of descriptors SDF data are provided
 * @author itetko
 *
 */

public class ChemPropServer extends SmilesOnlyAbstractServer
{
	/** Initializes supported task and workflow. */
	public ChemPropServer()
	{
		supportedTaskType = QSPRConstants.CHEMPROP;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);	
		batchApplySize = 10000;
		startPositionResults = 1;
	}

	@Override
	protected String getGPUCard(boolean forceCPU) { // always card 0, since environment is used to restrict only one card, which becomes 0
		return "" + (noGPU() || forceCPU ? NO_GPU : 0);
	}

	@Override
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean traningMode,
			boolean forceCPU) throws IOException, InterruptedException {
	}

	@Override
	String[] getCommands(ModelAbstractConfiguration receivedConf, boolean training) {
		ChemPropConfiguration cfg = (ChemPropConfiguration)receivedConf;

		String command [];

		if(training)
			command =  new String[] {executablePython, getExeFile(), 
					"--data_path", DATAFILE, 
					"--dataset_type", cfg.areClassificationData()?"classification":"regression",
							"--hidden_size", ""+ cfg.hidden,
							"--depth", ""+ cfg.depth,
							"--epochs", ""+ cfg.nepochs,
							"--batch_size", "" + cfg.batch,
							//"--graph_invariant_func", "rdkitunbranched-1-7-256",
							//"" + ((cfg.edges != null && cfg.edges)?"--virtual_edges":""),
							"--seed","" + cfg.getSeed(),
							"--split_sizes 0.9 0.1 0",
							"--save_dir",MODEL,
		};
		else
			command =  new String[] {executablePython, getAliasedPath() + "predict.py", 
					"--test_path", APPLYFILE, 
					"--checkpoint_dir",MODEL,
					"--preds_path",PREDICTIONS
		};

		if(!noGPU() && gpuCard != -1) {
			command = OCHEMUtils.append(command,"--gpu");
			command = OCHEMUtils.append(command,getGPUCard(false));
		}		
		return command;
	}


}
