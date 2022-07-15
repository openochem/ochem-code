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

import java.util.HashMap;

import com.eadmet.utils.FileUtils;

import qspr.metaserver.configurations.XGBOOSTConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LIBSVMUtils;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.ProcessUtils;
import qspr.workflow.utils.QSPRConstants;
import ml.dmlc.xgboost4j.java.Booster;
import ml.dmlc.xgboost4j.java.DMatrix;
import ml.dmlc.xgboost4j.java.XGBoost;

public class XGBOOSTServer extends MachineLearningExecutableAbstractServer {

	final String MODEL = "model"; 
	final String DATAFILE = "train.txt"; 
	final String APPLYFILE = "test.txt"; 

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors,
			LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf)
					throws Exception {

		XGBOOSTConfiguration configuration = (XGBOOSTConfiguration) receivedConf;

		LIBSVMUtils.writeLibSvmTrainingSet(dtDescriptors, dtExpValues, null, getAliasedFileName(DATAFILE));

		DMatrix trainMat = new DMatrix(getAliasedFileName(DATAFILE));

		HashMap<String, Object> params = new HashMap<String, Object>();
		if(configuration.lambda != null)params.put("lambda", configuration.lambda);
		params.put("eta", configuration.eta);
		params.put("max_depth", configuration.depth);
		params.put("silent", 1);
		params.put("objective",configuration.objective);
		params.put("seed", configuration.getSeed());
		if(ProcessUtils.getThreads()>1)params.put("nthread",ProcessUtils.getThreads());

		HashMap<String, DMatrix> watches = new HashMap<String, DMatrix>();

		//train a boost model
		Booster booster = XGBoost.train(trainMat, params, configuration.rounds, watches, null, null);
		booster.saveModel(getAliasedFileName(MODEL));

		configuration.storeModel(getAliasedFile(MODEL));

		return applyModel(dtDescriptors, configuration);
	}

	@Override
	protected DataTable applyModelBatch(DescriptorsTable dtDescriptors, ModelAbstractConfiguration receivedConf) throws Exception {

		FileUtils.saveBytesToFile(receivedConf.getSavedModelAsBytes(),getAliasedFileName(MODEL));

		LIBSVMUtils.writeLibSvmTestSet(dtDescriptors, getAliasedFileName(APPLYFILE));

		Booster booster = XGBoost.loadModel(getAliasedFileName(MODEL));
		DMatrix testMat = new DMatrix(getAliasedFileName(APPLYFILE));
		float[][] predicts = booster.predict(testMat);

		DataTable dtResults = new DataTable(true);
		dtResults.addColumn(QSPRConstants.PREDICTION_RESULT_COLUMN);
		dtResults.addColumn(QSPRConstants.DM+QSPRConstants.CLASSLAG);

		for(int i=0; i < dtDescriptors.getDataSize(); i++){
			AbstractDataRow row= dtResults.addRow();
			row.setValue(0, predicts[i][0]);
			if(receivedConf.areClassificationData())
				row.setValue(1, 0.5-Math.abs(0.5-predicts[i][0]));
		}

		return dtResults;
	}

	public XGBOOSTServer()
	{
		supportedTaskType = QSPRConstants.XGBOOST;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
		batchApplySize = 2*QSPRConstants.MODEL_REPOST_SIZE;
	}
}
