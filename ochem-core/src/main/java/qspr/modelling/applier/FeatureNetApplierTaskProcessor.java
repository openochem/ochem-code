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

package qspr.modelling.applier;

import java.io.Serializable;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.PendingTask.TaskType;
import qspr.exceptions.CalculationException;
import qspr.metaserver.configurations.ConsensusModelConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.DataReference;
import qspr.modelling.AbstractTaskProcessor;
import qspr.modelling.ModelApplierAttachment;
import qspr.modelling.ModelFactory;
import qspr.modelling.ModelProcessor;
import qspr.workflow.datatypes.WorkflowConfiguration;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

/**
 * A simplified applier TaskProcessor for FeatureNet and ConsensusModels
 * It reuses NoSQLReference object already created for these models to speed up analysis
 * @author itetko
 *
 */

public class FeatureNetApplierTaskProcessor extends AbstractTaskProcessor{
	DataReference data;
	Serializable configuration;
	public WorkflowNodeData wndResult;
	Model model;
	int dataSize;

	private static Logger logger = LogManager.getLogger(FeatureNetApplierTaskProcessor.class);

	public FeatureNetApplierTaskProcessor(DataReference data, Model model, PredictionScenario scenario, int molecules) throws Exception{
		this.data = data;
		dataSize = molecules;
		ModelProcessor processor = ModelFactory.getProcessor(model.template);
		processor.predictionScenario = scenario;
		processor.model = model;
		configuration = processor.getApplierConfiguration(true, false);
		if(configuration instanceof WorkflowConfiguration){
			WorkflowConfiguration conf = (WorkflowConfiguration) configuration;
			conf.repostSize = QSPRConstants.MODEL_REPOST_SIZE; // to fit to 2GB memory requirements			
		}
		this.model = model;
		setDescription = QSPRConstants.WEBSERVICE_TASK;
		taskName = "applying " + model.name;
		taskClass = TaskType.MODEL_APPLICATION;
	}

	@Override
	protected String getTaskType() {
		return configuration instanceof ConsensusModelConfiguration? QSPRConstants.CONSENSUS : QSPRConstants.Workflow;
	}

	@Override
	protected Serializable getTaskConfiguration() throws Exception {
		return configuration;
	}

	@Override
	protected Serializable getTaskData() throws Exception {
		return data;
	}

	/**
	 * Creating a fake attachment
	 * @throws Exception 
	 */

	@Override
	protected Serializable getAttachment() throws Exception {
		ModelApplierAttachment attachment = new ModelApplierAttachment();
		Basket basket = new Basket();

		for(int i=0 ; i < dataSize; i++)
			basket.entries.add(new BasketEntry(new ExperimentalProperty()));

		attachment.setWorkData(basket);
		return attachment;
	}

	@Override
	protected void onTaskReceived(Task task) throws Exception {
		this.taskReceived = true;
		if (task.isError())
		{

			logger.warn("Apply model task failed: "+task.getDetailedStatus());
			throw new CalculationException("Task failed: " + task.getDetailedStatus());
		}
		else if (task.isReady())
		{
			wndResult = (WorkflowNodeData)task.getResult();
			logger.info("Apply model has finished successfully");
		}
	}

	@Override
	protected void restoreFromAttachment(Serializable attachment) {
	}

	@Override
	protected Model getModel(){
		return model;
	}

}
