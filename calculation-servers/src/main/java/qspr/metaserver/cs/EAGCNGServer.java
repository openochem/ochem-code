
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

import com.eadmet.exceptions.CriticalException;

import qspr.metaserver.configurations.EAGCNGConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.workflow.utils.QSPRConstants;


/**
 *  Requires specific handling, since instead of descriptors SDF data are provided
 * @author itetko
 *
 */

public class EAGCNGServer extends SmilesOnlyAbstractServer
{
	int outputs;

	/** Initializes supported task and workflow. */
	public EAGCNGServer()
	{
		supportedTaskType = QSPRConstants.EAGCNG;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);	
		startPositionResults = 1;
	}

	@Override
	public boolean isCritical(String message) {
		if(message.contains("'NoneType' object has no attribute 'GetAtoms'"))
			throw new CriticalException("'NoneType' object has no attribute 'GetAtoms': " + message);
		return super.isCritical(message);
	}

	@Override
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean traningMode, boolean forceCPU) throws IOException{

		EAGCNGConfiguration conf = (EAGCNGConfiguration) configuration;

		conf.sanitize = true;

		//if(conf.nepochs > 500)conf.nepochs = 500; //TODO remove once finished analysis

		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		writer.write("[Task]");
		writer.write("\nmethod =  " + ("" + conf.method).toLowerCase());
		writer.write("\ntrain_mode = " + (traningMode?"True":"False"));
		writer.write("\nmodel_file = " + MODEL);
		writer.write("\ntask_names = ");
		for(int i = 0; i< conf.OutputValues(); i++)
			writer.write(QSPRConstants.PREDICTION_RESULT_COLUMN + i + ((i != (conf.OutputValues() -1))?",":""));
		writer.write("\ntrain_data_file = " + DATAFILE);
		writer.write("\napply_data_file = " + APPLYFILE);
		writer.write("\nresult_file = " + PREDICTIONS);
		if(conf.labelWeighting != null)
			writer.write("\nweights = " + WEIGHTS);

		writer.write("\n\n[Details]");
		writer.write("\nregression = " + (configuration.areClassificationData()? "False":"True"));
		writer.write("\ngpu = "+ getGPUCard(forceCPU));
		writer.write("\nseed = " + conf.getSeed() );
		writer.write("\nn_epochs = " + conf.nepochs);
		writer.write("\nbatch_size = " + conf.batchsize);
		writer.write("\ndropout = " + conf.dropout);
		writer.write("\nlearning_rate = " + conf.learningRate);
		writer.write("\nweightDecay = " + conf.weightDecay);
		writer.write("\nmethod = " + conf.method);
		writer.write("\nn_sgc1 = " + conf.n_sgc1);
		if(conf.n_sgc2 != null && conf.n_sgc2.trim().length() > 0) writer.write("\nn_sgc2 = " + conf.n_sgc2);
		if(conf.n_sgc3 != null && conf.n_sgc3.trim().length() > 0) writer.write("\nn_sgc3 = " + conf.n_sgc3);
		writer.write("\nn_den = " + conf.n_den);
		writer.write("\nnorm = " + conf.normalisation.toString().toLowerCase());
		writer.write("\ngate = " + (conf.gate?"True":"False"));
		writer.write("\nactivation = " + conf.getActivationFunction());

		writer.write("\n");
		writer.close();
	}


}
