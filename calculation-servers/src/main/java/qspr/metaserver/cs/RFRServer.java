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

import qspr.dao.Various;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.RFRConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.util.ExecutableRunner;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

/**
 * High level implementation of an interface to R server calculations
 * @author itetko
 *
 */

public class RFRServer extends MachineLearningExecutableAbstractServer
{

	final static String DATAFILE = "training.csv";
	final static String PREDICTIONS = "results.csv";
	final static String R = "r.r";
	final static String MODEL = "model.dat";
	String METHOD = "randomForest";

	public RFRServer()
	{
		supportedTaskType = QSPRConstants.RFR;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues,
			ModelAbstractConfiguration receivedConf) throws Exception {
		saveTRAIN((RFRConfiguration) receivedConf);
		saveDataCSV(dtDescriptors, dtExpValues, DATAFILE);
		runPythonWithConda(new String[]{ExecutableRunner.findExecutable("R"),"--slave","-f",getAliasedFileName(R)}, MODEL, null);
		receivedConf.storeModel(getAliasedFile(MODEL));
		return applyModel(dtDescriptors, receivedConf);
	}

	@Override
	protected DataTable applyModelBatch(DescriptorsTable dtDescriptors, ModelAbstractConfiguration receivedConf) throws Exception {
		saveAPPLY((RFRConfiguration) receivedConf);
		saveDataCSV(dtDescriptors, null, DATAFILE);
		saveModelToFile(receivedConf, getAliasedFileName(MODEL));
		runPythonWithConda(new String[]{ExecutableRunner.findExecutable("R"),"--slave","-f",getAliasedFileName(R)}, PREDICTIONS, null);
		return readResultValuesSingle(PREDICTIONS, receivedConf.optionsNumber, 1);
	}

	private void saveTRAIN(RFRConfiguration configuration) throws IOException {
		BufferedWriter f = getAliasedBufferedWriter(R);
		f.write("library(\"" + METHOD + "\")\n");
		f.write("data <- read.csv(file = \"" + DATAFILE + "\", header = TRUE)\n");

		f.write("set.seed(" + configuration.getSeed() + ")\n");
		f.write("model <-" + METHOD + "( " + QSPRConstants.PREDICTION_RESULT_COLUMN + 0 + " ~ ., data = data");
		if(configuration.numTrees != null)
			f.write(", ntree = " + configuration.numTrees);
		if(configuration.numFeatures != null && configuration.numFeatures !=0 )
			f.write(", mtry = " + configuration.numFeatures);
		if(configuration.maxDepth != null && configuration.maxDepth !=0 )
			f.write(", maxnodes = " + configuration.maxDepth);
		f.write(")\n");

		f.write("saveRDS(model,\""+ MODEL+ "\")\n");
		f.close();
	}

	private void saveAPPLY(RFRConfiguration configuration) throws IOException {
		BufferedWriter f = getAliasedBufferedWriter(R);
		f.write("library(\"" + METHOD + "\")\n");
		f.write("model <- readRDS(\"" + MODEL + "\")\n");
		f.write("data <- read.csv(file = \"" + DATAFILE + "\", header = TRUE)\n");
		f.write("pred <-predict(model,data)\n");
		f.write("write.csv(pred, file =\"" + PREDICTIONS + "\")\n");
		f.close();
	}

	protected void saveDataCSV(DescriptorsTable original, LabelsTable vals, String file) throws IOException, InterruptedException{

		boolean skip[]=null;
		boolean keepNames = true;
		boolean smiles = false;

		BufferedWriter writer  = getAliasedBufferedWriter(file);

		if(keepNames) {
			//if(smiles) writer.write(SMILES);

			for(int n=0; n < original.getColumnSize(); n++) {
				if(smiles || n > 0) writer.write(",");
				writer.write(DESC + n);
			}
			if(vals != null)
				/*				if(this instanceof MultiLearningAbstractServer) 
					for(int i=0;i<vals.getColumnSize();i++)
						writer.write("," + QSPRConstants.PREDICTION_RESULT_COLUMN + i);
				else
				 */					writer.write("," + QSPRConstants.PREDICTION_RESULT_COLUMN + 0); // only one; we do not expect here multiple columns
			writer.write("\n");
		}


		int size = vals == null? original.getDataSize() : vals.getDataSize();

		for(int i =0; i < size ; i++){
			if(skip != null && skip[i])continue;
			if(smiles)writer.write(Various.molecule.convertToFormat(original.getRawRow(i).getAttachment(QSPRConstants.SMILES_ATTACHMENT).toString(), QSPRConstants.SMILESNOSTEREO));
			if(original.getColumnSize()>0) { // do we have descriptors ?
				String descr[] = original.getScaledValuesString(i);
				for(int n=0; n < original.getColumnSize(); n++) {
					if(smiles || n > 0) writer.write(",");
					writer.write(descr[n]);
				}
			}

			if(vals!=null) {
				/*				if(this instanceof MultiLearningAbstractServer) 
					for(String val : vals.getScaledValuesString(i))
						if(NumericalTable.MISSED_VALUES.equals(val))writer.write(",");
						else
							writer.write("," + val);
				else
				 */					writer.append("," + vals.getOneValueOnlyString(i));
			}
			writer.write("\n");
		}
		writer.close();
	}


}
