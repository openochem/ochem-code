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
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Map;

import com.eadmet.exceptions.CriticalException;

import qspr.exceptions.CalculationException;
import qspr.interfaces.DataDrivenConfiguration;
import qspr.metaserver.configurations.BaggingConfiguration;
import qspr.metaserver.configurations.CompressedObject;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.metaserver.configurations.ValidationConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.transport.DataReference;
import qspr.metaserver.transport.DataReferenceException;
import qspr.metaserver.transport.DataReferenceFactory;
import qspr.metaserver.transport.DataReferenceFactory.DataReferencerType;
import qspr.metaserver.util.MachineLearningTask;
import qspr.workflow.utils.QSPRConstants;
import qspr.metaserver.util.TrainingSetSplitter;
import qspr.metaserver.util.aggregator.AbstractAggregator;
import qspr.metaserver.util.aggregator.AveragingAggregator;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

/**
 * The Bagging server
 * Originally developed by @author midnighter, then severely refactored by @author itetko
 * 
 * @midnighter
 * Concerns:
 * Usage of global variables, which makes 
 * (a) a stateful server 
 * (b) an unclear flow of data
 * 
 * It is a better practice to make servers stateless (seemingly to all other servers)
 *
 */

public class BaggingServer extends ValidationAbstractServer {
	private List<BitSet> savedValidationRowNums = null;
	private List<Map<Integer, Integer>> savedTrainingSetRepetions = null;
	private WorkflowNodeData savedTrainingSetData = null;
	private List<Object> models = null;
	private DataReference trainingSetReference;

	public BaggingServer() {
		supportedTaskType = QSPRConstants.BAGGING;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}

	protected void clean(){
		models = null; // everything is zeroed now
		savedValidationRowNums = null;
		savedTrainingSetRepetions = null;
		savedTrainingSetData = null;
	}

	@Override
	protected MachineLearningTask getValidationTask(WorkflowNodeData wdInput, ValidationConfiguration configuration, TrainingSetSplitter splitter, int bagNum)
			throws IOException {

		// logically -- there is no mapping, thus there is no training set and this is the application task
		if (splitter == null)
			return getApplyTask(wdInput, (BaggingConfiguration) configuration, bagNum); // in case if application

		MachineLearningTask vTask = new MachineLearningTask(this, wdInput, configuration, bagNum);

		Map<Integer, Integer> trainingSetMolecules = splitter.createBag(configuration.getBaggingType() == ValidationConfiguration.STRATIFIED,
				((BaggingConfiguration) configuration).getInstances());

		vTask.createTrainingAndValidationSets(trainingSetMolecules, wdInput, splitter, true);

		return vTask;
	}

	private MachineLearningTask getApplyTask(WorkflowNodeData wdInput, BaggingConfiguration configuration, int bagNum) throws IOException {

		MachineLearningTask vTask = new MachineLearningTask(this, wdInput, configuration, bagNum);

		ModelAbstractConfiguration conf = (ModelAbstractConfiguration) configuration.models.get(bagNum); // restore the model. It should be of this type!

		vTask.configuration = conf;

		if (vTask.configuration instanceof DataDrivenConfiguration) // if required, restore the training set too
		{
			Map<Integer, Integer> repeat = savedTrainingSetRepetions == null ? null : savedTrainingSetRepetions.get(bagNum);
			//CompactedWorkflowNodeData compactedDataset = new CompactedWorkflowNodeData(trainingSetReference, savedValidationRowNums.get(bagNum), repeat);
			//((DataDrivenConfiguration) vTask.configuration).setTrainingSet(compactedDataset);
			((DataDrivenConfiguration) vTask.configuration).setTrainingSet(extractDataset(trainingSetReference, savedValidationRowNums.get(bagNum), repeat));
		}
		vTask.wndInput = new WorkflowNodeData(wdInput.ports.get(0));
		vTask.keepModel = false;

		return vTask;
	}

	public WorkflowNodeData extractDataset(DataReference originalDatasetRef, BitSet excludedRows, Map<Integer, Integer> repetitionsMap) throws DataReferenceException {

		WorkflowNodeData originalDataset = (WorkflowNodeData) DataReferenceFactory.createReferencer(originalDatasetRef).getReference(originalDatasetRef);

		WorkflowNodeData extractedDataset = originalDataset.getEmptyCopy();
		int savedSetSize = originalDataset.ports.get(1).getRowsSize();

		for (int row = 0; row < savedSetSize; row++)
			if (!excludedRows.get(row))
			{
				int repetitions = 1;
				if (repetitionsMap != null && repetitionsMap.containsKey(row))
					repetitions = repetitionsMap.get(row);
				for (int i = 0; i < repetitions; i++) {
					for (int k = 0; k < originalDataset.ports.size(); k++)
						extractedDataset.ports.get(k).addRow(originalDataset.ports.get(k).getRow(row));
				}
			}

		return extractedDataset;
	}
	@Override
	protected DataTable trainModel(DescriptorsTable dtInput, LabelsTable dtLabels, ModelAbstractConfiguration receivedConfiguration) throws Exception {
		return trainModel(dtInput,dtLabels.getRawData(),(BaggingConfiguration) receivedConfiguration);
	}


	protected DataTable trainModel(DescriptorsTable dtInput, DataTable dtLabels, BaggingConfiguration configuration) throws Exception {

		DataTable dtResults = trainValidationModels(dtInput, dtLabels, configuration);

		if (savedValidationRowNums != null && savedValidationRowNums.size() > 0) {
			// Model is data-driven. Save the training set to the bagging configuration
			configuration.validationRowNumbers = new CompressedObject<List<BitSet>>();
			configuration.validationRowNumbers.set(savedValidationRowNums);
			configuration.trainingSet = new CompressedObject<WorkflowNodeData>();
			configuration.trainingSet.set(new WorkflowNodeData(dtInput.getSlice(0, dtLabels.getRowsSize()).getRawData()).addPort(dtLabels));
			configuration.trainingSetRepetions = new CompressedObject<List<Map<Integer, Integer>>>();
			configuration.trainingSetRepetions.set(savedTrainingSetRepetions);
		}
		configuration.models = models; // save model

		if (configuration.keepIndividualPredictions != null && configuration.keepIndividualPredictions)
			updateIndividualPredictions(dtInput,configuration,dtResults);

		clean(); // cleaning all temporal arrays
		return dtResults;
	}


	private void updateIndividualPredictions(DescriptorsTable dtInput, BaggingConfiguration configuration, DataTable dtResults) throws Exception{
		PredictionScenario savedScenario = null;

		if(configuration.models==null)return; // if models are not saved, individual predictions cannot be saved!

		if(configuration.models!=null)for(Object mod:configuration.models){
			if (mod instanceof ModelAbstractConfiguration) {
				ModelAbstractConfiguration m = ((ModelAbstractConfiguration) mod);
				if(savedScenario==null)savedScenario=m.predictionScenario;

				// we can store all scenarios and restore them for each model
				// However, according to our current code only one PredictionScenario is used in apply model
				if(savedScenario!=m.predictionScenario)
					throw new CriticalException("This code should be modified to work with different prediction scenarios");
				m.predictionScenario=PredictionScenario.DISTANCE_ONLY;
			}
		}

		DataTable predictionResults = applyModel(dtInput, configuration);
		for (int col = 0; col < dtResults.getColumnsSize(); col++)
			if (dtResults.getColumn(col).startsWith(QSPRConstants.INDIVIDUAL_PREDICTIONS))
				for (int row = 0; row < dtResults.getRowsSize(); row++)
					dtResults.setValue(row, col, predictionResults.getValue(row, col));

		for(Object mod:configuration.models ){
			if (mod instanceof ModelAbstractConfiguration)
				((ModelAbstractConfiguration) mod).predictionScenario=savedScenario;
		}
	}


	@Override
	protected DataTable applyModelBatch(DescriptorsTable dtDescriptors, ModelAbstractConfiguration receivedConf) throws Exception {

		BaggingConfiguration configuration = (BaggingConfiguration) receivedConf;

		if (configuration.trainingSetRepetions != null)
			savedTrainingSetRepetions = configuration.trainingSetRepetions.get(); // uncompressing only one time to speed-up processing of bagging tasks
		if (configuration.trainingSet != null)
		{
			savedTrainingSetData = configuration.trainingSet.get(); // uncompressing only one time to speed-up processing of bagging tasks
			trainingSetReference = DataReferenceFactory.createReferencer(DataReferencerType.MEMORY).saveReference(savedTrainingSetData,"memory");
		}

		if (configuration.validationRowNumbers != null)
			savedValidationRowNums = configuration.validationRowNumbers.get();

		return applyValidationModels(dtDescriptors, configuration);
	}

	@Override
	void saveModel(MachineLearningTask bagTask) throws IOException, ClassNotFoundException, CalculationException {

		if (bagTask.keepModel) {

			if (savedValidationRowNums == null)
				savedValidationRowNums = new ArrayList<BitSet>();
			if (savedTrainingSetRepetions == null)
				savedTrainingSetRepetions = new ArrayList<Map<Integer, Integer>>();
			if (models == null)
				models = new ArrayList<Object>();

			if (bagTask.getWndOutput().ports.get(1).getRowsSize() > 0) { // we may not have cfg if we are in doNotSave model mode
				Object modelConfiguration = bagTask.getWndOutput().ports.get(1).getValue(0, 0);
				if (modelConfiguration instanceof DataDrivenConfiguration) {
					((DataDrivenConfiguration) modelConfiguration).removeTrainingSet();
					savedValidationRowNums.add(listToBitSet(bagTask.validationRowNums));
					savedTrainingSetRepetions.add(bagTask.trainingRowRepetitions);
				}

				if (modelConfiguration instanceof ModelAbstractConfiguration)
					((ModelAbstractConfiguration) modelConfiguration).predictionScenario = null;

				models.add(modelConfiguration);
			}
		}
		bagTask.clean();
	}

	@Override
	AbstractAggregator getAggregator(int size, ValidationConfiguration conf) 
	{
		boolean keepPredictions = (conf.keepIndividualPredictions != null && conf.keepIndividualPredictions);
		return  new AveragingAggregator(size, keepPredictions);
	}

	private BitSet listToBitSet(List<Integer> list) {
		BitSet b = new BitSet();
		for (Integer i : list)
			b.set(i);
		return b;
	}
/*
	@Override
	protected int overRun(ValidationConfiguration configuration) {
		BaggingConfiguration conf = (BaggingConfiguration)configuration;
		return conf.noOverrun !=null && conf.noOverrun? 0:  (int)Math.ceil(conf.ensembleSize/10.);
	}
*/
}
