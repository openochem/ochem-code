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
import java.util.Map;

import qspr.exceptions.CalculationException;
import qspr.metaserver.configurations.CrossValidationConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.ValidationConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.util.MachineLearningTask;
import qspr.metaserver.util.TrainingSetSplitter;
import qspr.metaserver.util.aggregator.AbstractAggregator;
import qspr.metaserver.util.aggregator.AveragingAggregator;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

public class CrossValidationServer extends ValidationAbstractServer
{
	MachineLearningTask fullsetTask;

	public CrossValidationServer()
	{
		supportedTaskType = QSPRConstants.CROSSVALIDATION;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		setInputFlowGroup(1);
		setInputFlowGroup(2);
		setInputFlowGroup(3);
	}

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtLabelsOriginal, ModelAbstractConfiguration receivedConfiguration)
			throws Exception
	{
		CrossValidationConfiguration configuration = (CrossValidationConfiguration) receivedConfiguration;
		configuration.ensembleSize = configuration.ensembleSize + 1; // to calculate also task with all molecules

		DataTable dtResults = trainValidationModels(dtDescriptors, dtLabelsOriginal.getRawData(), configuration);

		if (configuration.saveModels())// only if we need it!
			configuration.taskConfiguration = (ModelAbstractConfiguration)fullsetTask.getWndOutput().ports.get(1).getValue(0, 0); // this is how we will send back configuration

		return dtResults;
	}

	@Override
	void saveModel(MachineLearningTask bagTask) throws IOException, ClassNotFoundException, CalculationException
	{
		// we are interested only in the last model!
		if (bagTask.keepModel)
			fullsetTask = bagTask;
		else
			bagTask.clean();
	}

	@Override
	AbstractAggregator getAggregator(int size, ValidationConfiguration conf)
	{
		return new AveragingAggregator(size, false);
	}

	@Override
	protected MachineLearningTask getValidationTask(WorkflowNodeData wdInput, ValidationConfiguration config, TrainingSetSplitter splitter, int taskNum)
			throws IOException
	{

		MachineLearningTask vTask = new MachineLearningTask(this, wdInput, config, taskNum);

		Boolean fullSetTask = taskNum == (config.ensembleSize - 1);

		if (fullSetTask && config.getBaggingType() != ValidationConfiguration.STRATIFIED)
			vTask.replicateTrainingSet(wdInput);
		else {
			Map<Integer, Integer> trainingSetMoleculeHashes = splitter.createFold(config.getBaggingType() == ValidationConfiguration.STRATIFIED, taskNum,
					config.ensembleSize - 1);
			vTask.createTrainingAndValidationSets(trainingSetMoleculeHashes, wdInput, splitter, fullSetTask);
		}

		vTask.keepModel = fullSetTask; // we do not need models from validation loops or if this is specifically requested

		return vTask;
	}

	@Override
	protected DataTable applyModelBatch(DescriptorsTable dtDescriptors, ModelAbstractConfiguration receivedConf) throws Exception
	{
		throw new IOException(" Cross-validation Task should not be used in apply mode!");
	}

}
