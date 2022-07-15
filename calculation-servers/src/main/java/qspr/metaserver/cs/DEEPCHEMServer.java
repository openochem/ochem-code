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
import java.io.IOException;

import com.eadmet.utils.OSType;

import qspr.metaserver.configurations.DeepChemConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.utils.QSPRConstants;


/**
 *  Requires specific handling, since instead of descriptors SDF data are provided
 * @author itetko
 *
 */

public class DEEPCHEMServer extends SmilesOnlyAbstractServer
{
	String methodtype = null;

	@Override
	protected CONDA getCondaEnvironment(){
		return OSType.isMac()?CONDA.RDKIT:null;
	}

	/** Initializes supported task and workflow. */
	public DEEPCHEMServer()
	{
		supportedTaskType = QSPRConstants.DEEPCHEM;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
		batchApplySize = QSPRConstants.MODEL_REPOST_SIZE/3 + 1;
	}

	@Override
	protected String methodType(){
		return methodtype;
	}

	
	@Override
	public int saveAggregatedData(String filename, DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws Exception {
		int records = super.saveAggregatedData(filename, dtDescriptors, dtExpValues, receivedConf);

		if(dtExpValues == null) return records;
		BufferedReader br = exeRunner.getAliasedBufferedReader(filename);

		int length=0;
		String line;
		while ((line = br.readLine()) != null) {
			if(line.contains(QSPRConstants.PREDICTION_RESULT_COLUMN)) continue;
			String tokens[] = line.split(",");
			if(tokens[0].length()>length)length=tokens[0].length();
		}
		DeepChemConfiguration conf = (DeepChemConfiguration) receivedConf;
		conf.maxlength = length;
		br.close();

		return records;
	}
	
	
	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean traningMode, boolean forceCPU) throws IOException{

		DeepChemConfiguration conf = (DeepChemConfiguration) configuration;

		conf.sanitize = true;

		methodtype = conf.isTEXTCNN()?"TEXTCNN"+(conf.maxlength == null?"":conf.maxlength):supportedTaskType.toUpperCase();

		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		writer.write("[Task]\n");
		writer.write("\ntrain_mode = " + (traningMode?"True":"False"));
		writer.write("\nmodel_file = " + MODEL);
		writer.write("\ntask_names = ");
		for(int i = 0; i< conf.OutputValues(); i++)
			writer.write(QSPRConstants.PREDICTION_RESULT_COLUMN + i + ((i != (conf.OutputValues() -1))?",":""));
		if(traningMode){
			writer.write("\ntrain_data_file = " + DATAFILE);
			if(conf.labelWeighting != null)
				writer.write("\nweights = " + WEIGHTS);
		}
		else
		{
			writer.write("\napply_data_file = " + APPLYFILE);
			writer.write("\nresult_file = " + PREDICTIONS);
		}


		writer.write("\n\n[Details]");

		String model = "";
		switch(conf.method) {
		case DAG:
			model = "DAG";
			break;
		case GRAPH_CONV:
			model = "GraphConv";
			break;
		case MPNN:
			model = "MPNN";
			break;
		case TEXTCNN:
			model = "TextCNN";
			break;
		case WEAVE:
			model = "Weave";
			break;
		default:
			break;
		}

		writer.write("\ndeepchem_model =  " + model);
		writer.write("\nseed = " + conf.getSeed() );
		writer.write("\ngpu = "+ getGPUCard(forceCPU));

		writer.write("\ntrain_proportion = " + (1-conf.early));
		writer.write("\nn_epochs = " + conf.nepochs);
		//writer.write("\nlearning_rate = " + conf.learning_rate);
		//writer.write("\ndropout = " + conf.dropout);
		//writer.write("\nn_embeddings = 10");

		//writer.write("\ndense_layer_size = " + conf.dense_layer_size);
		//writer.write("\ngraph_conv_layers = " + conf.graph_conv_layers);
		//writer.write("\nM = " + conf.M);
		//writer.write("\nT = " + conf.T);
		//writer.write("\nn_hidden = " + conf.n_hidden);
		//writer.write("\nn_embedding = " + conf.n_embedding);

		//ModelAbstractConfiguration.DataType datatype = conf.getExtendedDataType();
		//Boolean classtype = conf.isClassification();
		//writer.write("\nregression = " + 
		//		!((datatype == ModelAbstractConfiguration.DataType.CLASSIFICATION || datatype == ModelAbstractConfiguration.DataType.MULTICLASSIFICATION) && classtype)  
		//		);
		writer.write("\n");
		writer.close();
	}
	
	/*

	void saveConfig(ModelAbstractConfiguration configuration, DescriptorsTable mols, boolean traningMode, boolean forceCPU) throws IOException{

		DeepChemConfiguration conf = (DeepChemConfiguration) configuration;

		methodtype = conf.isTEXTCNN()?"TEXTCNN"+(conf.maxlength == null?"":conf.maxlength):supportedTaskType.toUpperCase();

		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		writer.write("[Task]\n");
		writer.write("\nseed = " + conf.getSeed() );
		writer.write("\ndeepchem_model =  " + ("" + conf.method).toLowerCase());
		writer.write("\ntrain_mode = " + traningMode);
		writer.write("\nmodel_save_dir = ./model");
		writer.write("\ntask_names = ");
		for(int i = 0; i< conf.OutputValues(); i++)
			writer.write(QSPRConstants.PREDICTION_RESULT_COLUMN + i + ((i != (conf.OutputValues() -1))?",":""));
		writer.write("\ntrain_data_file = " + DATAFILE);
		if(!traningMode) {
			writer.write("\ntest_data_file = " + APPLYFILE);
			writer.write("\ntest_result_file = " + PREDICTIONS);
		}else
			if(conf.labelWeighting != null)
				writer.write("\nweights = " + WEIGHTS);

		writer.write("\n\n[Details]");
		writer.write("\ngpu = "+ getGPUCard(forceCPU));

		writer.write("\ntrain_proportion = " + (1-conf.early));
		writer.write("\nnum_epochs = " + conf.nepochs);
		writer.write("\nlearning_rate = " + conf.learning_rate);
		writer.write("\ndropout = " + conf.dropout);
		writer.write("\ndense_layer_size = " + conf.dense_layer_size);
		writer.write("\ngraph_conv_layers = " + conf.graph_conv_layers);
		writer.write("\nM = " + conf.M);
		writer.write("\nT = " + conf.T);
		writer.write("\nn_hidden = " + conf.n_hidden);
		writer.write("\nn_embedding = " + conf.n_embedding);

		ModelAbstractConfiguration.DataType datatype = conf.getExtendedDataType();
		Boolean classtype = conf.isClassification();
		writer.write("\nregression = " + 
				!((datatype == ModelAbstractConfiguration.DataType.CLASSIFICATION || datatype == ModelAbstractConfiguration.DataType.MULTICLASSIFICATION) && classtype)  
				);
		writer.write("\n");
		writer.close();
	}

*/
}
