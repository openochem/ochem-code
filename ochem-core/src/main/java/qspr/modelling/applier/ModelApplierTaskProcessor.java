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

import java.io.IOException;
import java.io.Serializable;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Hibernate;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.PendingTask.TaskType;
import qspr.exceptions.CalculationException;
import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.workflow.utils.QSPRConstants;
import qspr.modelling.AbstractTaskProcessor;
import qspr.modelling.ApplicabilityDomain;
import qspr.modelling.ModelApplierCache;
import qspr.modelling.ModelFactory;
import qspr.modelling.ModelStatistics;
import qspr.modelling.SetStatistics;
import qspr.util.AccessChecker;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.exceptions.UserFriendlyException;

/**
 * The class responsible for application of a single model.
 * It is based on a concept of basket (can be also a virtual one), which contains molecules to be predicted
 * It is normally created and used only by ModelApplier.
 * 
 * @author midnighter
 *
 */
public class ModelApplierTaskProcessor extends AbstractTaskProcessor
{
	private ModelApplier parent;
	public Model model;

	public WorkflowNodeData wndResult;
	public DataTable dtDescriptors;
	public int propertiesCount = 0;

	private ModelApplierCache cache;

	/**
	 * Indicate whether to use cache for this task
	 */
	private Boolean useCache;

	private Map<String, ApplicabilityDomain> applicabilityDomains = new HashMap<String, ApplicabilityDomain>();


	public int getCachedCount(){
		return cache == null?0:cache.getCachedCount();
	}

	public ModelApplierTaskProcessor(ModelApplier parent, Model model)
	{
		try{
			model.checkIfCanBeApplied();
		}catch(Exception e){
			throw new UserFriendlyException(e.getMessage());
		}
		this.taskClass = TaskType.MODEL_APPLICATION;
		allowFailures = true;
		this.parent = parent;
		this.model = model;
		setDescription = parent.compoundsProvider.setDescription;
	}

	public void setUseCache(){
		useCache = true;
	}

	@Override
	public void start()
	{
		parent.compoundsProvider.status.addListener(statusTracker);

		if(setDescription == null)setDescription = parent.compoundsProvider.setDescription;

		if (model.id != null)
		{
			initialiseModel();
			propertiesCount = model.modelMappings.size();
			Hibernate.initialize(model.modelMappings);
			model = Globals.initializeAndUnproxy(model);
		}
		super.start();
	}

	public void addValidationSet() throws Exception
	{
		if (wndResult == null)
			return;

		boolean consistent = model.isStoredConsistently();

		parent.fetchRecordsFromDatabase();
		// Update statistics

		Basket basket = parent.compoundsProvider.getBasket();

		for (ModelMapping mmm : model.modelMappings)
		{
			ModelMapping mm = (ModelMapping) Globals.session().get(ModelMapping.class, mmm.id);

			ModelStatistics ms = ((ModelStatistics)mm.statisticsRecalculated.getObject());
			SetStatistics ss = new SetStatistics(basket, mm);
			wndResult.ports.get(0).reset();
			ss.setPredictions(wndResult.ports.get(0), mm, null);

			ss.setId = QSPRConstants.VALIDATION + model.getValidationSets().size();
			if (ms.getSetStatisticsByBasket(basket.id) != null)
			{
				if (ms.getSetStatisticsByBasket(basket.id).setId.startsWith(QSPRConstants.TRAINING))
					throw new UserFriendlyException("Cannot use the training set as a validation set");
				// We have already something for this set. Replace it with the new data.
				ss.setId = ms.getSetStatisticsByBasket(basket.id).setId;
				ms.sets.remove(ms.getSetStatisticsByBasket(basket.id));
			}
			ms.sets.add(ss);

			mm.statisticsRecalculated.updateObject();

			if (consistent)
				mm.statisticsOriginal.setObject(mm.statisticsRecalculated.getObject());

			Globals.session().saveOrUpdate(mm);
		}

		model = (Model) Globals.session().get(Model.class, model.id);
		model.addValidationSet(basket);
		Globals.session().saveOrUpdate(model);
	}

	@Override
	public void onTaskReceived(Task task) throws IOException, CalculationException
	{
		this.taskReceived = true;
		if (task == null)
		{
			// Task can be null, which means that we have had everything cached
			if (getCachedCount() > 0)
				wndResult = new WorkflowNodeData(cache.getCachedResult(ModelApplierCache.LOSTRESULT));
			else
				throw new UserFriendlyException("Prediction cache is empty - internal error");
		}
		else
		{
			if (task.isError())
			{
				// If we have something cached - we use it rather than failing completely.
				if (getCachedCount() > 0)
					wndResult = new WorkflowNodeData(cache.getCachedResult(task.getDetailedStatus()));
				else
				{
					logger.warn("Model "+model.name+" failed: "+task.getDetailedStatus());
					throw new CalculationException("Task failed: " + task.getDetailedStatus());
				}
			}
			else if (task.isReady())
			{
				long start = Calendar.getInstance().getTimeInMillis();

				wndResult = (WorkflowNodeData)task.getResult();
				initialiseModel();
				model.setColumnsNames(wndResult.ports.get(0).getColumns());

				if (!task.isError()){
					DataTable res = wndResult.ports.get(0);

					logger.info("Profiling: Getting task took " + (Calendar.getInstance().getTimeInMillis() - start)+"ms"); start = Calendar.getInstance().getTimeInMillis();

					if (cache != null){
						cache.cachePredictions(res);
						wndResult.ports.set(0,cache.mergeCachedResult(res)); // reusing received results

						logger.info("Profiling: Caching task results took " + (Calendar.getInstance().getTimeInMillis() - start)+"ms"); start = Calendar.getInstance().getTimeInMillis();

						task.setResult(wndResult); // updating all results with new ones

						logger.info("Profiling: Saving task  took " + (Calendar.getInstance().getTimeInMillis() - start)+"ms"); start = Calendar.getInstance().getTimeInMillis();

					}
				}
				logger.info("Model "+model.name+" has finished successfully");
			}
		}
	}

	public ApplicabilityDomain getApplicabilityDomain(String dmName, Integer mappingId) throws Exception
	{
		String hash = "" + dmName + "," + mappingId;
		if (applicabilityDomains.get(hash) != null)
			return applicabilityDomains.get(hash);

		if (wndResult == null || !ApplicabilityDomain.hasDM(model))
			return null;

		if (dmName == null)
			dmName = ModelStatistics.get(model.modelMappings.get(mappingId)).sets.get(0).distancesToModel.get(0);

		String dmColumnName = model.modelMappings.size() > 1 
				? "DM" + mappingId + ":" + dmName 
						: QSPRConstants.DM_PREFIX + dmName;	

		if (!wndResult.ports.get(0).containsColumn(dmColumnName))
			return null; // there is no needed DM in the predictioins table

		ApplicabilityDomain applicabilityDomain = new ApplicabilityDomain();
		applicabilityDomain.setModel(model.modelMappings.get(mappingId), dmName, null);
		applicabilityDomain.setTargetModelResults(this);
		applicabilityDomains.put(hash, applicabilityDomain);

		return applicabilityDomain;
	}

	public void resetApplicabilityDomains(String dmName) throws Exception
	{
		applicabilityDomains.clear();
		for (int i = 0; i < model.modelMappings.size(); i++)
			getApplicabilityDomain(dmName, i);		
	}

	private static Logger logger = LogManager.getLogger(ModelApplierTaskProcessor.class);

	@Override
	protected String getTaskType()
	{
		return QSPRConstants.Workflow;
	}

	@Override
	protected Serializable getTaskConfiguration() throws Exception
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected Serializable getTaskData() throws Exception
	{
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public Task createTask() throws Exception
	{
		logger.info("Starting creation of an applied task for " + model.name);
		Globals.restartAllTransactions(true);

		model = Model.unproxyModel(model);

		Basket basket = prepareBasket();

		if (basket.entries.size() > 0)
		{

			initialiseModel();

			setStatus("Preparing the calculation task for posting...");
			Task task = ModelFactory.createModelApplierTask(model, basket, true, parent.forceUpdateDescriptorCache, parent.scenario, parent.defaultConditions, parent.repostSize);
			task.debug = taskDebugLevel = parent.taskDebugLevel;

			logger.info("Applier task has been created "+model.name);
			return task;
		}

		if (cache != null)
		{
			setStatus("Checking already cached results ...");
			wndResult = new WorkflowNodeData(cache.getCachedResult(ModelApplierCache.LOSTRESULT));
			logger.info("All predictions for model " + model.name + " were cached");
			return null;
		}

		throw new UserFriendlyException("No molecules to predict were provided.");
	}

	private Basket prepareBasket(){

		if (useCache == null)
			useCache = parent.useCache && model.published;

		Basket basket = parent.compoundsProvider.getBasket();
		basket.evict(false);
		model = Model.unproxyModel(model);

		if (basket.entries.size() == 0)
			throw new UserFriendlyException("No molecules to predict were provided.");

		logger.info("Cached entries are analysed for " + model.name);

		if (useCache || parent.scenario != PredictionScenario.PREDICTION_ONLY) // only for published models or if CV are required
		{
			useCache = true;
			setStatus("Retrieval of previously calculated and cached predictions...");
			cache = new ModelApplierCache(model, basket, parent.scenario, parent.defaultConditions);
			return cache.getNonCachedEntriesBasket();
		}

		logger.info("Basket is ready for " + model.name);

		return basket;
	}

	@Override
	public int getTaskPriority()
	{
		int priority = parent.defaultTaskPriority != null ? parent.defaultTaskPriority : TaskPriority.HIGH;
		return priority = Math.min(priority, AccessChecker.getMaximumTaskPriority(Globals.getCurrentUser()));
	}

	@Override
	protected Serializable getAttachment()
	{
		// Do not store a heavyweight attachment until we got a saved pending task
		if (pTask == null)
			return null;
		return parent.getAttachment().setCache(cache);
	}

	@Override
	protected Model getModel()
	{
		return model;
	}

	@Override
	protected String getSetDescription()
	{
		return setDescription;
	}

	@Override
	protected void restoreFromAttachment(Serializable _attachment)
	{
	}

	public void setCache(ModelApplierCache cache) {
		this.cache = cache;		
	}

	@Override
	public int bonusMultiplier() {
		return QSPRConstants.APPLIER_BONUS;
	}

	public void initialiseModel() {
		model = Model.initialiseModel(model);
	}
}

