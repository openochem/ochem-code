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

import qspr.metaserver.configurations.HAMNETConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.workflow.utils.QSPRConstants;

public class HAMNETServer extends SmilesOnlyAbstractServer{
	/** Initializes supported task and workflow. */
	public HAMNETServer()
	{
		supportedTaskType = QSPRConstants.HAMNET;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);	
		batchApplySize = 10000;
		startPositionResults = 1;
		MODEL = "model.pt";
	}

	@Override
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean train, boolean forceCPU) throws IOException{

		HAMNETConfiguration conf = (HAMNETConfiguration) configuration;

		conf.sanitize = true;
		
		File f = new File(getAliasedFileName(CFG));
		if(f.exists()) {
			File ff = new File(getAliasedFileName(CFG+".train"));
			if(!ff.exists())
				f.renameTo(ff);
		}

		BufferedWriter writer = getAliasedBufferedWriter(CFG);

		writer.write("[Task]\n"+
				"train_mode = "+  (train?"True":"False")+ "\n"+
				(train?"\ntrain_data_file = " + DATAFILE:"apply_data_file = " + APPLYFILE) + "\n"+
				"result_file = " + PREDICTIONS + "\n" +
				"model_file = "+ MODEL + "\n");

		writer.write("task_names = ");

		for(int i = 0; i< conf.OutputValues(); i++)
			writer.write(QSPRConstants.PREDICTION_RESULT_COLUMN + i + ((i != (conf.OutputValues() -1))?",":""));				

		writer.write("\n\n[Details]\n"
				+   "batch_size = " + conf.batch_size  + "\n"
				+ 	"gpu = "+ getGPUCard(forceCPU)   + "\n");

		if(train)writer.write(
				"seed = "+ conf.getSeed() + "\n"
						+   "n_epochs = "  + conf.nepochs  + "\n"
						+   "learning_rate = "  + conf.learningRate  + "\n"
						+   "dropout = 0.3\n"
						+	"train_proportion = " + (1. - conf.early) +"\n"  
						+   "n_best_nets = 10\n"
						+   "restart = False\n"
						+ (conf.isRefined()?"refine = True":"") + "\n"
				);

		writer.write("\n");
		writer.close();
		return;
	}

}
