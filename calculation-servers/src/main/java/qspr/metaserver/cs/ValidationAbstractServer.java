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
import java.util.Calendar;

import qspr.exceptions.CalculationException;
import qspr.metaserver.configurations.CrossValidationConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.ValidationConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.util.MachineLearningTask;
import qspr.metaserver.util.TrainingSetSplitter;
import qspr.metaserver.util.aggregator.AbstractAggregator;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.CalculationTask;
import qspr.workflow.utils.CalculationTaskSet;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.CriticalException;

/**
 * 
 * @author itetko
 * Abstract class to work with validation protocols
 * 
 */

public abstract class ValidationAbstractServer extends MachineLearningAbstractServer
{

	protected DataTable trainValidationModels(DescriptorsTable dtMolecules, DataTable dtLabelsOriginal, ValidationConfiguration configuration)
			throws Exception
	{
		setStatus("Running Validation server, tasks received used memory: " + Integer.toString((int) (usedMemory() / (1024 * 1024))) + " MB");

		if(checkValuesUniqueness(dtLabelsOriginal))throw new CriticalException("All property (target) values are the same and equal "+ 
				dtLabelsOriginal.getValue(0, 0) + ". The model can't be built: you can't develop a model for molecules with identical (or without) property values.");

		if(dtLabelsOriginal==null)throw new IOException(" dtLabelsOriginal==null -- use applyValidationModel instead");

		CalculationTaskSet taskSet = new CalculationTaskSet(this);

		// We should use taskConfiguration and not Validation Configuration
		// Only task Configuration contains information about classes of compounds
		ModelAbstractConfiguration conf = (ModelAbstractConfiguration) configuration.taskConfiguration;
		LabelsTable taskLabels = new LabelsTable(dtLabelsOriginal, conf); // this class is used to define type of labels for data splitting

		TrainingSetSplitter splitter = 	new TrainingSetSplitter(dtMolecules, taskLabels, out, configuration.mixtureValidation, conf.seed==null?QSPRConstants.SEED:conf.seed);
				
		if (conf.getExtendedDataType() != ModelAbstractConfiguration.DataType.CLASSIFICATION && conf.getExtendedDataType() != ModelAbstractConfiguration.DataType.MULTICLASSIFICATION  
				&& configuration.getBaggingType() == ValidationConfiguration.STRATIFIED)
			throw new CriticalException("Stratified Validation can be only applied to classification task with two classes.");

		WorkflowNodeData wndInput = new WorkflowNodeData(dtMolecules.getRawData());
		wndInput.addPort(dtLabelsOriginal);

		ModelAbstractConfiguration taskConfiguration = (ModelAbstractConfiguration) configuration.taskConfiguration;

		// Because of ambiguity, information whether we need to save models can be in several places
		taskConfiguration.saveModels = configuration.saveModels = configuration.saveModels() && taskConfiguration.saveModels() ? null : false;

		int ensemble = configuration.ensembleSize + overRun(configuration);

		for (int i = 0; i < ensemble; i++)
		{
			long timeStart = Calendar.getInstance().getTimeInMillis();

			MachineLearningTask bagTask = getValidationTask(wndInput, configuration, splitter, i);

			postTask(taskSet, bagTask, i == ensemble - 1 || allowLocalCalculations); // last task will be calculated on the server itself

			out.println("Task has been posted, it took " + (Calendar.getInstance().getTimeInMillis() - timeStart) / 1000 + " sec");
			setStatus("Running Validation server, task " + i + " posted  (out of " + ensemble + ") used memory: "
					+ Integer.toString((int) (usedMemory() / (1024 * 1024))) + " MB");
		}

		setStatus("Running Validation server, tasks posted used memory: " + Integer.toString((int) (usedMemory() / (1024 * 1024))) + " MB");

		taskSet.calculate(true,overRun(configuration));

		setStatus("Running validation server, ready to aggregate final results; used memory: " + Integer.toString((int) (usedMemory() / (1024 * 1024))) + " MB");

		return aggreateResults(dtLabelsOriginal.getRowsSize(), dtMolecules.getDataSize(), taskSet, configuration);

	}


	/**
	 * Specifies number of additional bags to run in order to skip them and finish calculations before waiting for long servers
	 * @param ensemble
	 * @return
	 */
	protected int overRun(ValidationConfiguration configuration){
		return 0;
	}


	/**
	 * Check whether all values in the table are the same
	 */

	static boolean checkValuesUniqueness(DataTable data) {

		if(!data.isCompactRowFormat()) return false; // no check is done for non-compact format

		if(data.getRowsSize() < 2) return true;

		for(int col=0; col<data.getColumnsSize();col++){
			double val0 = (Double) data.getValue(0, col);
			for(int i =1; i< data.getRowsSize(); i++){
				double val1 = (Double) data.getValue(i, col);
				if(val1 != val0) return false;
			}
		}
		return true;
	}

	protected DataTable applyValidationModels(DescriptorsTable dtMolecules, ValidationConfiguration configuration)
			throws Exception
	{
		setStatus("Applying Validation server, tasks received used memory: " + Integer.toString((int) (usedMemory() / (1024 * 1024))) + " MB");

		CalculationTaskSet taskSet = new CalculationTaskSet(this);

		WorkflowNodeData wndInput = new WorkflowNodeData(dtMolecules.getRawData());

		if (dtMolecules.getDataSize() > 0)
		{
			for (int i = 0; i < configuration.ensembleSize; i++)
			{
				long timeStart = Calendar.getInstance().getTimeInMillis();

				MachineLearningTask bagTask = getValidationTask(wndInput, configuration, null, i);

				postTask(taskSet, bagTask, true); // only local calculations will be tried first

				out.println("Task " + bagTask.taskId+ " has been posted, it took " + (Calendar.getInstance().getTimeInMillis() - timeStart) / 1000 + " sec");

				if(bagTask.taskId<0)  // otherwise the task was posted to metaserver (e.g., we do not have enough memory or local server is not availabel)
					while(!bagTask.isReady())
						Thread.sleep(100);

				setStatus("Applying Validation server, task " + i + " posted  (out of " + configuration.ensembleSize + ") used memory: "
						+ Integer.toString((int) (usedMemory() / (1024 * 1024))) + " MB");
			}

			setStatus("Applying Validation server, tasks posted used memory: " + Integer.toString((int) (usedMemory() / (1024 * 1024))) + " MB");

			taskSet.calculate(true);
		}

		return aggreateResults(0, dtMolecules.getDataSize(), taskSet, configuration);

	}


	private DataTable aggreateResults(int trainingSetSize, int wholeSet, CalculationTaskSet taskSet, ValidationConfiguration conf) throws Exception
	{
		AbstractAggregator aggregator = getAggregator(wholeSet, conf);

		int count = 0, bag = 0;

		// Gather statistics for all the points from different bagging tasks
		for (CalculationTask task : taskSet.tasks)
		{
			out.println("Aggregating results: bag  " + bag);

			MachineLearningTask bagTask = (MachineLearningTask) task;

			if(bagTask.isReady() && bag < conf.ensembleSize){

				DataTable predictions = bagTask.getWndOutput().ports.get(0);

				int bagTrainSize = bagTask.trainingSetSize;
				int bagValidationInternal = bagTask.validationRowNums.size();

				out.println("Aggregating results: bag  " + bag + " initial trainingSetSize=" + trainingSetSize + " bagTrainingSetSize=" + bagTrainSize
						+ " internalBagSize=" + bagValidationInternal + " external bag=" + (predictions.getRowsSize() - bagTrainSize - bagValidationInternal));

				// first we aggregate rows for internal validation set
				for (int k = bagTrainSize; k < bagTrainSize + bagValidationInternal; k++)
				{
					aggregator.aggregateRow(bagTask.validationRowNums.get(k - bagTrainSize), predictions, k, bag, conf.ensembleSize);
					count++;
				}

				out.println("aggregated " + count + " internal validation results");

				// now we aggregate rows (if any) for the prediction set
				for (int k = bagTrainSize + bagValidationInternal; k < predictions.getRowsSize(); k++)
				{
					aggregator.aggregateRow(k - bagTrainSize - bagValidationInternal + trainingSetSize, predictions, k, bag, conf.ensembleSize);
					count++;
				}
				out.println("aggregated " + count + " all validation results");

				saveModel(bagTask);
				bag++;
			}else{
				out.println("skipping model ready:"+task.isReady()); // skipping over-run models
				task.kill(); // killing task to avoid any further calculations
			}
		}

		out.println("Aggregated " + count + " rows in total");

		try {
			return aggregator.getAggregatedResult(conf instanceof CrossValidationConfiguration);
		}catch(Exception e) {
			throw new CriticalException("Calculation failed due to insufficient number of samples in validation folds." + 
					(conf instanceof CrossValidationConfiguration ? " Decrease number of folds " + 
							(conf.getBaggingType() == ValidationConfiguration.STRATIFIED ? " or/and disable the stratified validation." :"") :" Try to increase the number of folds."));
		}

	}

	abstract AbstractAggregator getAggregator(int size, ValidationConfiguration conf);

	abstract void saveModel(MachineLearningTask bagTask) throws IOException, ClassNotFoundException, CalculationException;

	abstract protected MachineLearningTask getValidationTask(WorkflowNodeData wdInput, ValidationConfiguration config, TrainingSetSplitter mappings, int i)
			throws IOException;


}
