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

import java.io.Serializable;
import java.util.HashSet;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.configurations.ConsensusModelConfiguration;
import qspr.metaserver.configurations.ConsensusModelConfiguration.IndividualModel;
import qspr.metaserver.configurations.DescriptorType;
import qspr.metaserver.configurations.DescriptorsApplyModelConfiguration;
import qspr.metaserver.cs.util.WebServiceCalculationTask;
import qspr.workflow.utils.QSPRConstants;
import qspr.metaserver.util.aggregator.AveragingAggregator;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.CalculationTask;
import qspr.workflow.utils.CalculationTaskSet;

/** A server to train and apply consensus models based on a number of different individual models / by Midnighter
 * So far, consensus models are just averaging aggregators
 */

public class ConsensusModelServer extends WorkflowNodeServer
{
	protected static final Logger logger = LogManager.getLogger(ConsensusModelServer.class);

	@Override
	public WorkflowNodeData calculate(WorkflowNodeData task, Serializable configuration) throws Exception
	{
		DataTable dtMolecules = task.ports.get(0);

		if(dtMolecules.getRowsNoErrorsSize() == 0)
			return new WorkflowNodeData(new DataTable());

		ConsensusModelConfiguration config = (ConsensusModelConfiguration) configuration;

		CalculationTaskSet taskSet = calculateTaskSet(task, this, config);

		setStatus(" have been finished.");

		int properties = 0, models = taskSet.tasks.size();
		DataTable predictions[] = new DataTable[models];
		for (int m = 0 ; m < taskSet.tasks.size(); m++)
		{
			CalculationTask cTask = taskSet.tasks.get(m);
			DataTable dtResult = ((WebServiceCalculationTask)cTask).getWndOutput().ports.get(0);
			predictions[m] = dtResult;
			properties = dtResult.getColumnsSize() > properties ? dtResult.getColumnsSize() : properties;
			if(dtResult.containsMissedDescriptors()) config.allowErrors = true;
		}

		setStatus("Calculations have been received: aggregating results.");

		properties = properties/2; // since we store both values and accuracies

		AveragingAggregator aggregator = new AveragingAggregator(1, config.type, config.allowErrors, true); // just for one row
		aggregator.stdDmName = QSPRConstants.CONSENSUS_STD ; // to differentiate from BAGGING-STD

		float predictedValues[] = new float[models]; // by number of models
		float predictedErrors[] = new float[models];

		DataTable dtResult = null;

		List<HashSet<Long>> modelsForProperties = config.getModelsForProperties();

		for (int row = 0; row < dtMolecules.getRowsSize();row++){ // by number of molecules

			aggregator.cleanRows(); // each row - one property
			String error = null;

			for(int property = 0 ; property < properties; property++){

				for (int model = 0 ; model < models; model++){ // by all models

					DataTable dtRes = predictions[model];

					if(dtRes.getRow(row).isError()){
						predictedValues[model] = predictedErrors[model] = Float.NaN;
						if(error == null)error = "Predictions failed for model #" + config.individualModels.get(model).id + " " + dtRes.getRow(row).detailedStatus;
						continue;
					}

					predictedValues[model] = ((Double)dtRes.getValue(row, property * 2)).floatValue();

					if(modelsForProperties != null && !modelsForProperties.get(property).contains(config.individualModels.get(model).id))
						predictedValues[model] = Float.NEGATIVE_INFINITY;					

					switch(config.type){
					case OPTIMAL:
					case AVERAGE: break; // accuracy is not required!
					case BEST_MODEL:
					case WEIGHTED_AVERAGE:
						predictedErrors[model] = ((Double)dtRes.getValue(row, property * 2 + 1)).floatValue(); 
						if(config.units.get(property).equals(QSPRConstants.CLASS)) predictedErrors[model] = 1 - predictedErrors[model]; // using error and not accuracy
						break;
					case RMSE_WEIGHTED:
						predictedErrors[model] = config.weights.get(model)[property]; // using saved values
					}
				}

				if(config.classes != null && config.classes.get(property)>2) // fooling the Aggregator
					AveragingAggregator.substituteWithMaxClass(predictedValues, predictedErrors, config.classes.get(property));

				aggregator.addPredictions(0, property, predictedValues, predictedErrors);
			}

			DataTable tab = aggregator.getAggregatedResult(false);
			AbstractDataRow oneRow = tab.getRow(0);
			if(oneRow.isError() && error != null)
				oneRow.setError(error);

			if(dtResult == null) dtResult = tab;
			else
				dtResult.addRow(oneRow);

		}

		return new WorkflowNodeData (dtResult);
	}

	static public CalculationTaskSet calculateTaskSet(WorkflowNodeData wndInput, WorkflowNodeServer server, ConsensusModelConfiguration consConf) throws Exception
	{
		CalculationTaskSet taskSet = new CalculationTaskSet(server);

		for (IndividualModel iModel : consConf.individualModels)
		{
			DescriptorsApplyModelConfiguration config = new DescriptorsApplyModelConfiguration(iModel.id);
			config.scenario = consConf.predictionScenario;
			config.units = consConf.units.<String>toArray(new String[]{});

			DescriptorType deskType = new DescriptorType(config);

			deskType.skipCache = true; // we do NOT cache Consensus predictions -- otherwise too many and actually there is no sense...

			WebServiceCalculationTask cTask = new WebServiceCalculationTask(deskType, wndInput, server, false);

			taskSet.addTask(cTask);
		}

		taskSet.post();
		taskSet.calculate(false);

		return taskSet;
	}

	public ConsensusModelServer()
	{
		supportedTaskType = QSPRConstants.CONSENSUS;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}

}

// A wrapper to apply models through OCHEM web-services
