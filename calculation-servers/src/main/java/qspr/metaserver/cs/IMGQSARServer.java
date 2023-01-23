
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

import com.eadmet.exceptions.CriticalException;
import com.eadmet.utils.OSType;

import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.IMGQSARConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.util.AliasedDirectory;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.utils.QSPRConstants;

/**
 * Requires specific handling, since instead of descriptors SDF data are
 * provided
 * 
 * @author itetko
 *
 */

public class IMGQSARServer extends SmilesOnlyAbstractServer {
	/** Initializes supported task and workflow. */
	public IMGQSARServer() {
		supportedTaskType = QSPRConstants.IMGQSAR;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
		batchApplySize = 10000;
		startPositionResults = 1;
	}

	protected CONDA getCondaEnvironment(){
		return CONDA.RDKIT;
	}

	@Override
	public boolean isCritical(String message) {
		if(message.contains("Explicit valence for atom"))
			throw new CriticalException("Explicit valence for atom: " + message);
		return super.isCritical(message);
	}
	
	@Override
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean train, boolean forceCPU) throws IOException{

		File f = getAliasedFile("");
		File model = new File(OSType.isMac()?"/etc/ochem/cddd":"/etc/source/imgqsar/");
		
		AliasedDirectory.alias(model, model.list(), f);
		
		IMGQSARConfiguration conf = (IMGQSARConfiguration) configuration;

		conf.sanitize = true;

		f = new File(getAliasedFileName(CFG));
		if(f.exists()) {
			File ff = new File(getAliasedFileName(CFG+".train"));
			if(!ff.exists())
				f.renameTo(ff);
		}

		String tasks = "task_names = ";
		for(int i = 0; i< conf.OutputValues(); i++)
			tasks += QSPRConstants.PREDICTION_RESULT_COLUMN + i + ((i != (conf.OutputValues() -1))?",":"");				

		int augx = conf.augmentation == null || conf.augmentation == 1? 0: conf.augmentation;
		int augy = conf.augmentApplySet == null? 1 : conf.augmentApplySet;
		
		if(augx  > 0) augx = -2;
		
		String configur=
				"[Task]\n"
						+ "train_mode = " + (train? "True":"False") +"\n"
						+ "model_file = " + MODEL+"\n"
						+ (train?"train_data_file = " + DATAFILE : "apply_data_file = "+ APPLYFILE) + "\n"
						+ tasks +"\n"
						+ "result_file= " + PREDICTIONS + "\n"
						+ "\n"
						+ "[Details]\n"
						+ "batch_size = " + conf.batch_size +"\n"
						+ "gpu = " + getGPUCard(forceCPU) + "\n"
						+ "n_epochs = " + conf.nepochs+ "\n"
						+ "learning_rate = "+ conf.learning_rate +"\n" //**
						+ "augment= " + augx +" \n"
						+ "apply_augment = " + augy +" \n"
						+ "n_best_nets = 5\n" //**
						+ "train_proportion= " + (1 - conf.early) +" \n"
						+ "seed = " + conf.getSeed() + "\n" 
						//+ (conf.backbone != null?"backbone = " + conf.backbone + "\n":"")
						//+ "shuffle = " + (conf.shuffleAllData()?"True":"False") + "\n"
						;

		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		writer.write(configur );
		writer.close();
		return;
	}

}
