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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.JAXBException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.type.LongType;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.Globals;
import qspr.business.ModelPeer;
import qspr.business.PendingTaskPeer;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.Property;
import qspr.entities.Unit;
import qspr.export.ExportableModel;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.util.WrapperThread;
import qspr.workflow.utils.QSPRConstants;

/**
 * Starts multiple models given a number of model configuration templates.
 * This is one of the core classes for the Comprehensive Modeling utility
 * 
 * @author midnighter
 *
 */
public class MultipleModelsStarter extends WrapperThread
{
	private static transient final Logger logger = LogManager.getLogger(MultipleModelsStarter.class);

	public Long trainingSetId;
	public List<Long> validationSetId;

	/**
	 * The task priority for the calculation tasks
	 */
	public int taskPriority = TaskPriority.LOW;

	/**
	 * Model configuration templates
	 */
	public List<ExportableModel> eModels = new ArrayList<ExportableModel>();

	public int skipped = 0, started = 0, noncompatible = 0;
	public boolean skipDublicates = true;
	//public boolean cancel;
	public Map<Property, Unit> units = new HashMap<Property, Unit>();

	public void wrapped() throws InterruptedException, JAXBException
	{
		setStatus("Initialising...");

		long entries = Repository.basket.countEntries(trainingSetId);

		try{
			if(entries == 0)
				throw new UserFriendlyException("Training set does not have any records. Models cannot be calculated.");
			Repository.user.checkEligibility(entries, QSPRConstants.MODEL_BONUS);
		}catch(UserFriendlyException e){
			setStatus("Error: " + e.getMessage());
			return;
		}

		try
		{
			startAllModels();
			String status = "Finished: ";
				status += (started - skipped - noncompatible) + " models have been started";
			if (skipped > 0)
				status += ", " + (skipped) + " already existing models";
			if(noncompatible > 0) status += ", and " + (noncompatible) + " non-compatible models";

			if(skipped + noncompatible > 0)
				status +=  " have been skipped";

			setStatus(status + ".");
		}
		catch (Exception e)
		{
			if(e.getMessage() !=null && !e.getMessage().startsWith("Error:")){
				e.printStackTrace();
				setStatus("Error: " + e.getMessage());
			}else
				setStatus(e.getMessage());
		}
	}

	private void startAllModels() throws Exception
	{
		started = skipped = noncompatible = 0;

		int MAX_MODELS = QSPRConstants.MAX_MODELS;

		Basket b = (Basket) Globals.session().get(Basket.class, trainingSetId);

		PendingTaskPeer.updateTaskStatuses(Globals.userSession().user, b, 200);

		int tasks = PendingTaskPeer.countTasks(Globals.userSession().user, b, Task.ASSIGNED);
		tasks += PendingTaskPeer.countTasks(Globals.userSession().user, b , Task.INIT);

		if(tasks > MAX_MODELS * 4)
			throw new UserFriendlyException("There are currently " + tasks + " models, which are calculating for \"" + b.name + "\" set. Wait untill these calculations are finished before starting new ones.");

		int all = 0;

		for (ExportableModel eModel : eModels){

			setStatus("Starting model " + (++started) + " out of " + eModels.size());

			if(((CDSConfiguration) eModel.attachment.configuration).hasDescriptors() && all >= MAX_MODELS) { 
				skipped++;
				continue;
			}

			eModel.trainingSetId = trainingSetId;
			if (validationSetId != null && !validationSetId.isEmpty())
				eModel.validationSetId = validationSetId;

			Model model = eModel.createModel();

			if (cancelRequested)
				return;

			model.configurationHash = eModel.getMD5();

			if (skipDublicates && alreadyExists(model))
			{
				logger.info("Skipping a model, it already exists");
				skipped++;
				continue;
			}

			if(!model.isCompatibleModelAndDescriptors()) {
				logger.info("Skipping a model, it cannot start in such combination of method/descriptors");
				noncompatible++;
				continue;
			}

			all++;

			model.session = Globals.userSession();
			model.createMappings(units);
			model.name = ModelPeer.getModelName(model, false);
			if (cancelRequested)
				return;
			ModelProcessor processor = ModelFactory.getProcessor(model.template);
			processor.model = model;
			processor.model.defaultTaskPriority = taskPriority;
			processor.updateSessionTime = false;

			Globals.commitAllTransactions(); // do not keep the transaction open
			processor.start();
			processor.exitAfterPosting();
			processor.join();
			if(QSPRConstants.ERROR_STATUS.equals(processor.model.status))
				throw new UserFriendlyException(processor.model.detailedStatus !=null && processor.model.detailedStatus.startsWith("Error:") ? processor.model.detailedStatus : "Error: "+processor.model.detailedStatus);
			Globals.startAllTransactions();
		}
	}


	private boolean alreadyExists(Model model)
	{
		model.attachment.updateObject();
		Object a = Globals.session().createSQLQuery("select count(*) c from Model where training_set_id=:tsId and configuration_hash=:hash")
				.addScalar("c", LongType.INSTANCE)
				.setParameter("tsId", model.trainingSet.id)
				.setParameter("hash", model.configurationHash)
				.uniqueResult();

		return ((Long) a) > 0;
	}
}
