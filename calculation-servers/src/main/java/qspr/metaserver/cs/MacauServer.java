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
import java.io.Writer;

import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.OSType;

import qspr.metaserver.configurations.MACAUConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

public class MacauServer extends MultiLearningAbstractServer{

	private static final String CFG = "config.ini";
	final String MODEL =  "macau.tar"; 
	final String PREDICTIONS = "results.csv";

	final String applyModel = "apply.py";
	final String trainModel = "train.py";

	final String DESCRIPTORS = "descrs.csv";
	final String PROPERTIES = "prop.csv";

	boolean saveValues;

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors,
			LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf)
					throws Exception {
		DescriptorsTable dtOriginal = dtDescriptors;

		MACAUConfiguration conf = (MACAUConfiguration) receivedConf;

		saveConfig(conf);

		saveAggregatedData(getAliasedFileName("/"),dtDescriptors, dtExpValues, receivedConf);
		String [] commands = {"/opt/conda/envs/rxnfp/bin/python",getExeFile(),"--training"}; // N.B.! FIX works only for linux for very specific setting
		executeBinaryBash(commands, MODEL);

		receivedConf.storeModel(getAliasedFile(MODEL));

		return applyModel(dtOriginal, receivedConf);
	}

	private void saveConfig(MACAUConfiguration conf) throws IOException {
		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		writer.write("[macau]\n");
		writer.write("\ntest_percent = " + conf.test_percent );
		writer.write("\nnum_latent  =  " + conf.num_latent);
		if(!conf.containAlsoRegressionData())writer.write("\nprecision = probit");
		else
			if(conf.isAdaptive())writer.write("\nprecision = adaptive");
			else
				writer.write("\nprecision = " + 1./conf.accuracy);
		writer.write("\nburnin = " + conf.burnin);
		writer.write("\nsamples  = " + conf.samples);
		writer.write("\n");
		writer.close();
	}

	@Override
	protected DataTable applyModelMulti(DescriptorsTable dtDescriptors, MultiLearningAbstractConfiguration receivedConf) throws Exception {

		saveAggregatedData(getAliasedFileName("/"),dtDescriptors,null, receivedConf);

		saveModelToFile(receivedConf, getAliasedFileName(MODEL));

		String [] commands = {"/opt/conda/envs/rxnfp/bin/python",getExeFile(),"--prognosis"};

		runPythonWithConda(commands, PREDICTIONS,CONDA.RDKIT);
		DataTable res = readResultValues(PREDICTIONS,receivedConf);

		return res;
	}

	protected String[] addOptional(String[] commands,
			String additionalParam) {

		if(additionalParam != null)
			for(String par : additionalParam.split("\\s+"))
				commands = OCHEMUtils.append(commands,par);

		return commands;
	}


	@Override
	public int saveAggregatedData(String filename, DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws IOException {

		BufferedWriter bw = getAliasedBufferedWriter(filename + DESCRIPTORS);
		saveValues = false; // saving compressed descriptors
		for(int i=0; i<dtDescriptors.getColumnSize();i++)
			bw.write((i>0?",":"")+QSPRConstants.DESCRIPTOR+i);
		bw.write("\n");
		int savedRecords = saveAggregatedDescriptorsAndValues(dtDescriptors, dtExpValues, bw, receivedConf); // how many we can compress
		bw.close();

		if(dtExpValues!= null) { // saving also target values in another file
			saveValues = true;
			bw = getAliasedBufferedWriter(filename + PROPERTIES);
			for(int i=0; i<dtExpValues.getColumnSize();i++)
				bw.write((i>0?",":"")+QSPRConstants.PREDICTION_RESULT_COLUMN+i);
			bw.write("\n");
			saveAggregatedDescriptorsAndValues(dtDescriptors, dtExpValues, bw, receivedConf); // how many we can compress
			bw.close();
		}

		return savedRecords;

	}

	@Override
	protected void saveOneRow(int mol, DescriptorsTable dtDescriptors, String[] values, Writer writer) throws IOException {

		if(saveValues)
			for(int i=0; i <values.length; i++)
				writer.append((values[i] == null ? "": values[i]) + (i !=values.length -1? ",":""));
		else{
			String descs[]=dtDescriptors.getScaledValuesString(mol);
			for(int i = 0;i<descs.length;i++)
				writer.append(descs[i]+(i+1 !=descs.length?",":""));
		}

		writer.append(OSType.endLine());
	}

	public MacauServer()
	{
		supportedTaskType = QSPRConstants.MACAU;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}

}