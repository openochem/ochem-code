
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
import qspr.metaserver.configurations.TRANSNNConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.workflow.utils.QSPRConstants;


/**
 *  Requires specific handling, since instead of descriptors SDF data are provided
 * @author itetko
 *
 */

public class TRANSNNServer extends SmilesOnlyAbstractServer
{
	/** Initializes supported task and workflow. */
	public TRANSNNServer()
	{
		supportedTaskType = QSPRConstants.TRANSNN;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);	
		batchApplySize = 10000;
		startPositionResults = 0;
	}

	@Override
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean train, boolean forceCPU) throws IOException{

		TRANSNNConfiguration conf = (TRANSNNConfiguration) configuration;

		File f = new File(getAliasedFileName(CFG));
		if(f.exists()) {
			File ff = new File(getAliasedFileName(CFG+".train"));
			if(!ff.exists())
				f.renameTo(ff);
		}

		BufferedWriter writer = getAliasedBufferedWriter(CFG);

		writer.write("[Task]"+
				"\nmodel_file = "+ MODEL + "\n"
				);

		if(train) {
			writer.write("\ntrain_mode = True" +
					"\ntrain_data_file = " + DATAFILE +
					"\n\n[Details]\n" +
					"first-line=True\n" +
					"n_epochs = " + conf.nepochs + "\n"+
					"batch_size = " + conf.batch +"\n" +
					"early-stopping = " + (conf.early == 0?0:(1 - conf.early)) + "\n" +
					"learning_rate = " + conf.learning_rate + "\n" +
					"chirality = " + (conf.requiresStereochemistry()?"True":"False") + "\n" +
					"retrain = " + (conf.retrain== null || conf.retrain?"True":"False") + "\n" +
					(conf.fixedrate != null && conf.fixedrate? "fixed-learning-rate = True\n":"")
					);
		}
		else
			writer.write("train_mode = False\n"+
					"first-line=True\n" +
					"apply_data_file = " + APPLYFILE +"\n" +
					"result_file = " + PREDICTIONS + "\n" +
					"\n\n[Details]\n"
					);

		writer.write(
				"canonize = False\n" + 
						"gpu = "+ getGPUCard(forceCPU)  + "\n" + 
						(conf.random?"random = True\n":"")+
						"seed = " +conf.getSeed() +"\n"
				);

		writer.close();
		return;
	}

}
