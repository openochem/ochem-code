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

package qspr.modelling;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.entities.Basket;
import qspr.entities.ConditionSet;
import qspr.entities.Model;
import qspr.entities.ModelTemplate;
import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.metaserver.protocol.Task;
import qspr.workflow.utils.QSPRConstants;

public class ModelFactory 
{
	private static transient final Logger logger = LogManager.getLogger(ModelFactory.class);

	public static ModelProcessor getProcessor(String name)
	{
		if (name.equals(QSPRConstants.CONSENSUS))
			return new ConsensusModelProcessor();
		else if (name.equals(QSPRConstants.UPLOADED_MODEL))
			return new ModelUploadProcessor();		
		else{
			logger.info("Using default CDSModelProcessor for " + name);
			return new CDSModelProcessor(); // assuming this type of processor
		}
		//			throw new RuntimeException("Cannot work with model type " + name);
	}

	public static ModelProcessor getProcessor(ModelTemplate template)
	{
		return ModelFactory.getProcessor(template.name);
	}

	public static ModelProcessor getProcessor(Model model)
	{
		ModelProcessor processor = ModelFactory.getProcessor(model.template);
		processor.model = model;

		return processor;
	}

	public static Task createModelApplierTask(Model model, Basket basket, boolean recalculated, boolean forceRecalculateDescriptors, PredictionScenario scenario, 
			ConditionSet defaultConditions, int repostSize) throws Exception
	{
		// Generate applier configuration on the fly
		ModelProcessor processor = getProcessor(model);
		if (processor instanceof CDSModelProcessor)
			((CDSModelProcessor) processor).applierRepostSize = repostSize;
		processor.predictionScenario = scenario;
		processor.model = model;
		Task task = processor.createApplierTask();
		task.taskName = "Applying " + model.name;

		task.setConfiguration(processor.getApplierConfiguration(recalculated,forceRecalculateDescriptors));
		task.disableTaskLevelCache = true;
		if (basket != null)
			task.setData(processor.getDataForApplier(basket, defaultConditions));

		return task;
	}


}
