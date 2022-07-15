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

import java.io.IOException;
import java.io.Serializable;
import java.security.NoSuchAlgorithmException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Pattern;

import javax.xml.bind.JAXBException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.FlushMode;
import org.hibernate.Hibernate;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.business.ModelPeer;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.entities.Attachment;
import qspr.entities.AttachmentSource;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.PendingTask;
import qspr.entities.PendingTask.TaskType;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.metaserver.protocol.Mergable;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.DebugLevel;
import qspr.metaserver.util.MixtureAttachment;
import qspr.workflow.utils.QSPRConstants;
import qspr.util.BasicRecordMapper;
import qspr.util.WrapperThread;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.MAILERConstants;

abstract public class ModelProcessor extends WrapperThread
{
	public Model model;
	public BasicRecordMapper mapper;
	protected MyCalculationClient client = new MyCalculationClient();
	protected long time = Calendar.getInstance().getTimeInMillis();
	private DataTable dtPredictedValues;
	public PredictionScenario predictionScenario = PredictionScenario.PREDICTION_ONLY;
	public Integer defaultTaskPriority;
	public String preferredServer = OCHEMConfiguration.defaultPreferredServer;
	public int taskDebugLevel; // To post a task with the debug flag

	private boolean exitAfterPosting = false;

	public abstract Serializable getDataForTeacher() throws Exception;
	abstract Serializable getDataForApplier(Basket basket, ConditionSet defaultConditions) throws Exception;
	public abstract Serializable getApplierConfiguration(boolean recalculatedModel, boolean forceDescriptorRecalculation) throws Exception;
	abstract Serializable getTeacherConfiguration() throws Exception;
	abstract protected String getTaskType();

	static Integer post = 0;

	@Override
	public void wrapped() throws JAXBException 
	{

		time = Calendar.getInstance().getTimeInMillis();

		setStatus("Starting initialization ...");

		if (model.id != null)
			model = (Model) Globals.session().get(Model.class, model.id); // due to a strange behavior, merge -> get
		else
			Globals.session().saveOrUpdate(model);

		Repository.user.checkEligibility(Repository.basket.countEntries(model.trainingSet.id), QSPRConstants.MODEL_BONUS);

		model.attachment.updateObject();
		model.updateDescription();
		for (ModelMapping mm : model.modelMappings)
			Globals.session().saveOrUpdate(mm);

		String implicitValues = model.implicitValues;

		if(model.template.isDescriptorCalculationOnly())throw new UserFriendlyException("Descriptor calculation task is not allowed in this content. Please, report the bug to "+MAILERConstants.EMAIL_OCHEM);

		if (pTask == null)
			registerPendingTask(new PendingTask(TaskType.MODEL_TRAINING, 0)
					.setModel(model).setPriority(defaultTaskPriority != null ? defaultTaskPriority : model.defaultTaskPriority));
		else
		{
			pTask.taskId = 0;
			registerPendingTask(pTask);
		}
		model.taskId = 0;

		pTask.setPriority(TaskPriorityManager.getNewTaskPriority(pTask.getPriority(), pTask, model.session));

		String user = pTask.session.user.login;

		Globals.commitAllTransactions();

		time = Calendar.getInstance().getTimeInMillis();

		setStatus("Waiting in the queue for " + user + " ...");
		synchronized (post)
		{
			try {
				setStatus("Starting posting for " + user + " ...");
				Globals.startAllTransactions(); // Start the transaction AFTER we got into the synchronized part. That will prevent unnecessary timeouts. // Midnighter on Aug 22, 2011
				if (model.id != null)
					model = (Model) Globals.session().get(Model.class, model.id);

				model.implicitValues = implicitValues;

				prepare();

				if(!model.isCompatibleModelAndDescriptors())
					throw new UserFriendlyException("This model is incompatible with the configuration and cannot be started.");
				createEmptyStatistics();
				postTeacherTask();
				setStatus("Finished posting " + user + " ...");;
				model.taskId = pTask.taskId;  //TODO VERIFY THAT CODE is NOT BROKEN
				if (exitAfterPosting)
					return;

			}catch(Throwable e) {
				setStatus("Exception for " + user+ " " + e.getMessage());
				pTask.detailedStatus = e.getMessage();
				model.detailedStatus = e.getMessage();
				Globals.session().saveOrUpdate(pTask);
				return;
			}finally {
				Globals.commitAllTransactions();
			}
		}

		// Task is in standalone mode
		try
		{
			Task resultTask = waitForTask();
			setStatus("Finished calculations for " + user);;

			// We are back and task is finished, open a new transaction
			Globals.startAllTransactions();

			onTaskReceived(resultTask);
			saveModel();
			Globals.commitAllTransactions();
			setStatus(QSPRConstants.FINISHED);
		} catch (Throwable e)
		{
			exception = new IOException(e.getMessage());
			e.printStackTrace();
			model.status = QSPRConstants.ERROR_STATUS;
			model.detailedStatus = e.getMessage();
			setStatus("Error: " + e.getMessage());
			if (model.name != null) {
				if(!Globals.areTransactionsRunning())
					Globals.startAllTransactions();
				Globals.session().saveOrUpdate(model);
				Globals.commitAllTransactions();
			}
		}
	}

	public Task waitForTask() throws Exception
	{
		client.originalStatus = getStatus();
		client.setTolerateMetaserverDown();
		model.lastModification = new Timestamp(Calendar.getInstance().getTimeInMillis());
		Task task = null;
		while (task == null)
		{
			task = client.getTask(model.taskId);
			Thread.sleep(500);
		}
		model.lastModification = new Timestamp(Calendar.getInstance().getTimeInMillis());

		return task;
	}

	public void saveModel() throws Exception
	{
		//update model
		model.attachment.updateObject();
		if (model.readyModelAttachment != null)
		{
			model.readyModelAttachment.updateObject();
			logger.info("Saving the model, size " + model.readyModelAttachment.getDataLength() / (1024 * 1024) + "Mb");
		}
		model.microattachment.updateObject();
		model.updateDescription();
		Globals.session().saveOrUpdate(model);
		if (!model.template.isDescriptorCalculationOnly() && model.getModelData(true) == null)
			throw new UserFriendlyException("Warning: The model was not properly saved!");

		model = (Model) Globals.session().merge(model);

		updateStatistics();

		model.isStatisticsCalculated = true;
		model.setColumnsNames(dtPredictedValues.getColumns());
		Globals.session().saveOrUpdate(model);

		// disabling caching of all predictions: simply too many !
		//		model.clearCachedPredictions();
		//		model.saveCachedPredictions(dtPredictedValues);
		logger.info("Model " + model.name + " is saved");
	}

	private void updateStatistics() throws Exception {
		// Update statistics
		List<Set<Long>> moleculesInSets = new ArrayList<Set<Long>>();

		logger.info("updateStatistics: collectingImplicit values");

		boolean virtualFound = false;
		dtPredictedValues.reset();
		do {
			dtPredictedValues.nextRow();
			if(dtPredictedValues.getCurrentRow().getAttachment(QSPRConstants.IMPLICIT_ATTACHMENT) != null)virtualFound = true;
		}while(dtPredictedValues.hasMoreRows() && !virtualFound);

		if (virtualFound) 
			for (ModelMapping mm : model.modelMappings)
			{
				ModelStatistics msOriginal = (ModelStatistics) mm.statisticsOriginal.getObject();
				ModelStatistics msRecalculated = (ModelStatistics) mm.statisticsRecalculated.getObject();
				ModelStatistics msAnalysed = msOriginal.containsPredictions ? msRecalculated : msOriginal;

				for (int i = 0 ; i < msAnalysed.sets.size() ; i++ ){
					if(moleculesInSets.size() <= i)moleculesInSets.add(new HashSet<Long>());
					Set<Long> set = moleculesInSets.get(i);
					for(PointStatistics ps: msAnalysed.sets.get(i).points) 
						set.add(ps.id);
				}
				if(msAnalysed.sets.size() != moleculesInSets.size()) throw new Exception("Difefrent number of sets in ModelPammping");
			}

		logger.info("updateStatistics: collectingImplicit values finished");

		for (ModelMapping mm : model.modelMappings)
		{
			ModelStatistics msOriginal = (ModelStatistics) mm.statisticsOriginal.getObject();
			ModelStatistics msRecalculated = (ModelStatistics) mm.statisticsRecalculated.getObject();

			if (!model.template.isDescriptorCalculationOnly())
				if (!msOriginal.containsPredictions)
				{
					if (msRecalculated != null)
					{
						msOriginal = msRecalculated;
						mm.statisticsOriginal.setObject(msOriginal);
					}
					msOriginal.setPredictions(dtPredictedValues, mm, moleculesInSets);
					mm.statisticsOriginal.updateObject();
					mm.statisticsRecalculated.setObject(msOriginal);
				}
				else
					msRecalculated.setPredictions(dtPredictedValues, mm, moleculesInSets);

			mm.statisticsRecalculated.updateObject();
			mm.updateDescription();
			mm.model = model;
			Globals.session().saveOrUpdate(mm);
		}

		logger.info("updateStatistics: finished");

	}

	public void onTaskReceived(Task task) throws Exception
	{
		setStatus("Task received - processing the result");
		task.check();
		onTeacherFinished((WorkflowNodeData) task.getResult());
		model.attachment.updateObject();
		model.microattachment.updateObject();
		model.timeToComplete = (task.timeCompleted.getTime() - task.timeAssigned.getTime()) / 1000;
		model.recalculateModelSize();
	}

	protected void prepare() throws Exception
	{
		setStatus("Preparing to run");
		model.isStatisticsCalculated = false;

		setStatus("Initializing the training set...");
		model.initTrainingAndTestSetEntries();
		setStatus("Training set initialized");

		mapper = new BasicRecordMapper(model.trainingSet);
	}

	protected void createEmptyStatistics() throws NoSuchAlgorithmException
	{
		setStatus("Creating the empty statistics template");
		for (ModelMapping mm : model.modelMappings)
		{
			ModelStatistics ms = ModelStatistics.getEmptyStatistics(mm);
			model.recalculation = (mm.statisticsOriginal != null && ((ModelStatistics) mm.statisticsOriginal.getObject()).containsPredictions);
			if (mm.statisticsOriginal == null)
				mm.statisticsOriginal = new Attachment<ModelStatistics>(ms, AttachmentSource.ModelStatistics);
			mm.statisticsRecalculated = new Attachment<ModelStatistics>(ms, AttachmentSource.ModelStatistics);
		}
		setStatus("Empty statistics template created.");
	}

	protected void postTeacherTask() throws Exception
	{
		Serializable teacherConfiguration = getTeacherConfiguration();
		Task task = null;

		Globals.session().setFlushMode(FlushMode.COMMIT);

		if (model.taskId == null || model.taskId == 0)
		{
			task = createTeacherTask();
			if (teacherConfiguration instanceof Mergable)
				task.mergeConfiguration((Mergable) teacherConfiguration);
			else
				task.setConfiguration(teacherConfiguration);
			task.taskName = model.name;
			Serializable teacherData;

			preloadLazyData();
			Globals.commitAllTransactions();
			task.setData(teacherData = getDataForTeacher());
			Globals.restartAllTransactions(true);
			logger.info("Teacher data: " + teacherData);
			task.setPriority(pTask.getPriority());

			model.lastAccess = model.dateCreated = new Timestamp(Calendar.getInstance().getTimeInMillis());
			setStatus("Posting the task in teacher");

			task.debug = taskDebugLevel;

			if (preferredServer != null) {
				task.setPreferredServer(preferredServer);
				if(task.debug == DebugLevel.NONE)
					task.debug = DebugLevel.ALL;
			}

			client.setDeepSleepTime(10000); //TODO origin of this step is no clear; it is possibly related to bug with submission of multiple application tasks, which I have closed
			pTask.taskId = model.taskId = client.postTask(task);
			Globals.session().saveOrUpdate(pTask);
			logger.info("Saving pending task " + pTask.taskId + " " + pTask.model.name); // for debugging, remove when cleared

			setStatus("Saving model in the database");

			if (model.name == null)
				model.name = ModelPeer.getModelName(model, true);
			model.attachment.updateObject();

			model.microattachment.updateObject();

			setStatus("Model training");
			model.lastModification = new Timestamp(Calendar.getInstance().getTimeInMillis());

			Globals.session().saveOrUpdate(model);
			setStatus("Model training started");
		}

		for (ModelMapping modelM : model.modelMappings)
			Globals.session().saveOrUpdate(modelM);

		// Close transaction before going to a long-waiting-loop, otherwise connection timeouts are possible / Midnighter
		Globals.commitAllTransactions();

		// Free the memory....
		Basket basket = model.getFilteredSet(QSPRConstants.TRAINING);
		for (BasketEntry entry : basket.entries)
		{
			entry.ep.molecule.setData("");
			entry.ep.molecule.molImage = null;
		}
		System.gc();
		logger.info("Memory after GC: " + MemoryUtils.memorySummary());
	}

	protected void onTeacherFinished(WorkflowNodeData teacherResponse) throws Exception
	{
		dtPredictedValues = teacherResponse.ports.get(0);
	}

	protected void onApplierFinished(WorkflowNodeData applierResponse) throws Exception
	{
	}

	public static DataTable basketToSDFTable(Basket basket) throws Exception{
		return basketToSDFTable(basket, null);
	}

	/**
	 * Create a DataTable with molecules based on a basket.
	 * This function requires some refactoring. A more strict interface should be established for access/creation of the molecules datatables
	 */
	public static DataTable basketToSDFTable(Basket basket, DataTable dtInitial) throws Exception
	{
		Pattern pattern = Pattern.compile("END.+", Pattern.DOTALL);
		DataTable datatable;
		if (dtInitial == null)
		{
			datatable = new DataTable();
			datatable.id = "mols";
			datatable.addColumn(QSPRConstants.SDF_COLUMN);
		}
		else
			datatable = dtInitial;

		datatable.compressStrings = true;

		datatable.setColumnAttachment(QSPRConstants.SDF_COLUMN, "session-guid", Globals.userSession() == null ? null: Globals.userSession().guid);

		List<BasketEntry> entries = basket.entries;

		String message = QSPRConstants.ERROR_SMILES;

		for (BasketEntry entry : entries)
		{
			datatable.addRow();

			AbstractDataRow row = datatable.getCurrentRow();

			if (entry.ep.molecule == null)
			{
				row.setError(entry.ep.errorComment != null ? entry.ep.errorComment : message);
				continue;
			}

			// This is an in-memory molecule, reload it from DB
			if (entry.ep.molecule.mapping2 == null)
				entry.ep.molecule = Repository.molecule.getMolecule(entry.ep.molecule.id);

			if (entry.ep.molecule == null)
				row.setError("Empty molecule provided");
			else
				if (entry.ep.molecule.error != null)
					row.setError("Error in molecule: " + entry.ep.molecule.error);
				else
					if (!entry.ep.molecule.hasData())
						row.setError("Empty molecule provided");

			if (row.isError())
				continue;

			if (!entry.ep.molecule.isEmptyMolecule()) {
				String sdf = entry.ep.molecule.getData();
				try {
					int atoms = Various.molecule.getAtomCount(sdf);
					if(atoms >= QSPRConstants.MAX_ATOMS)
						row.setError("Molecule has " + atoms + " heavy atoms. OCHEM does not process molecules with >= " + QSPRConstants.MAX_ATOMS + " heavy atoms. ");
					else {
						sdf = pattern.matcher(sdf).replaceAll("END").replace("$$$$", "");
						datatable.setValue(QSPRConstants.SDF_COLUMN, sdf);
						row.addAttachment(QSPRConstants.INCHIKEYS, entry.ep.molecule.mapping1.inchi1);
					}
				}catch(Exception e) {
					row.setError("Failed to calculate atoms number or InChi - molecule cannot be processed");
				}
			}
			else
				datatable.setValue(QSPRConstants.SDF_COLUMN, null);

			if (entry.ep.molecule.mapping2.inchi2.length() == 32) // InChies are not available; molecule can be error
			{
				if (entry.ep.externalId == null)
					row.setError(entry.ep.errorComment != null ? entry.ep.errorComment : message);
			}
			/****/
			row.addAttachment(QSPRConstants.RECORD_ID_ATTACHMENT, entry.ep.id);

			if (entry.ep.externalId != null)
				row.addAttachment(QSPRConstants.EXTERNAL_ID, entry.ep.externalId);

			if (!entry.ep.molecule.isEmptyMolecule())
			{
				row.addAttachment(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM, entry.ep.molecule.mapping1.id.intValue());
				row.addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, entry.ep.molecule.mapping2.id.intValue());
				row.addAttachment(QSPRConstants.MOLECULE_ID_OCHEM, entry.ep.molecule.id);
			}
			else
			{
				if (entry.ep.externalId != null)
				{

					try
					{
						row.addAttachment(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM, Integer.valueOf(entry.ep.externalId.replace("UCB", "-")));
						row.addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, Integer.valueOf(entry.ep.externalId.replace("UCB", "-")));
					} catch (Exception e)
					{
						row.addAttachment(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM, entry.ep.id.intValue());
						row.addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, entry.ep.id.intValue());
					}

				}
				else
					row.setError("Empty molecule provided!");
			}

			try
			{
				addConditionMetadata(entry, row);
			}
			catch (Exception e)
			{
				e.printStackTrace();
				row.setError(e.getMessage());
			}

			if (datatable.currentRow % 10000 == 0 && datatable.currentRow > 0)
			{
				logger.info("Processed " + datatable.currentRow + " rows");
				if (MemoryUtils.getCurrentMemoryUsedFraction() > 0.95)
					throw new UserFriendlyException("Unfortunately, the dataset you're trying to model is too large and can not be handled now. Try using a smaller set or try again later.");
			}
		}

		if(datatable.getRowsSize() < basket.entries.size()) // should be at least the same number as in the basket can be more, since table can accumulate more entries
			throw new UserFriendlyException("Different number of molecules in the initial basket " + basket.entries.size() + " and in the SDF file with entries");

		return datatable;
	}


	//TODO to move to correct place
	private static void addConditionMetadata(BasketEntry entry, AbstractDataRow row) throws Exception
	{
		if (entry.ep.conditions == null) return; 

		boolean activeTx = Globals.session().getTransaction().isActive(); //FIXME: Remove all transaction-dependent code from here... somehow

		if (!activeTx)
			Globals.startAllTransactions();

		// Mixture Related Code
		ConditionSet conditions = (ConditionSet) Globals.session().get(ConditionSet.class, entry.ep.conditions.id);

		if (conditions.getValue(QSPRConstants.MIXTURE_CONDITION) != null)
		{
			MixtureAttachment ma = ExperimentalProperty.createMixtureAttachment(conditions.getValue(QSPRConstants.MIXTURE_CONDITION).textualValue);
			if(ma.fractions != null) 
				row.addAttachment(QSPRConstants.MIXTURE_ATTACHMENT, ma);
		}

		if (!activeTx)
			Globals.commitAllTransactions();

	}

	/**
	 * Preload all the lazy data (e.g., condition values).
	 * This is required to be able to close the connection safely before a long data preparation process
	 */
	protected void preloadLazyData()
	{
		logger.info("Initializing lazy data");
		Map<String, Basket> map = getAllSets();
		for (String setId : map.keySet())
		{
			for (BasketEntry be : model.getFilteredSet(setId).entries)
			{
				if (be.ep.conditions != null)
					Hibernate.initialize(be.ep.conditions.values);
			};
		}
	}

	/**
	 * Enumerate all the participating sets (training, validation, excluded)
	 */
	private Map<String, Basket> getAllSets()
	{
		Map<String, Basket> map = new HashMap<String, Basket>();
		Basket trainingSet = model.getFilteredSet(QSPRConstants.TRAINING);
		map.put(QSPRConstants.TRAINING, trainingSet);

		for (int i = 0; i < model.getValidationSets().size(); i++)
		{
			Basket validationSet = model.getFilteredSet(QSPRConstants.VALIDATION + i);
			if (validationSet != null)
				map.put(QSPRConstants.VALIDATION + i, validationSet);
		}

		Basket excludedSet = model.getFilteredSet(QSPRConstants.EXCLUDED);
		if (excludedSet != null)
			map.put(QSPRConstants.EXCLUDED, excludedSet);

		return map;

	}

	protected static DataTable basketToPropertiesTable(Basket basket, DataTable dtInitial, Model model, BasicRecordMapper mapper) throws Exception
	{
		DataTable dtExpValues;
		if (dtInitial == null)
		{
			dtExpValues = new DataTable(true);
			dtExpValues.id = "exp-values";
			dtExpValues.addColumn(QSPRConstants.VALUE);
			if (mapper.getRowsSize() > 1)
				dtExpValues.addColumn(QSPRConstants.CLASS);
		}
		else
			dtExpValues = dtInitial;

		List<BasketEntry> entries = basket.entries;
		List<ModelMapping> modelMappings = model.modelMappings;
		boolean handleRanges = model.attachment.getObject().datahandling.handleRanges;

		Double values[] = new Double[2];

		for (BasketEntry entry : entries)
		{
			dtExpValues.addRow();
			//converted value 

			try {
				for (ModelMapping modelMapping : modelMappings)
					if (modelMapping.matches(entry.ep))
					{
						if (modelMapping.property.isQualitative())
						{
							dtExpValues.setValue(new Double(model.getMappedOption(entry.ep.option.id)));
							if (dtExpValues.getValue() == null){
								dtExpValues.getCurrentRow().setError("No quantitative mapping specified for option " + entry.ep.option.name);
								continue;
							}
						}
						else
						{
							if (handleRanges)
							{
								dtExpValues.setValue(QSPRConstants.VALUE, values[0] = values[1] = entry.ep.getConvertedValue(modelMapping.unit));
								dtExpValues.getCurrentRow().addAttachment(QSPRConstants.PREDICATE_ATTACHMENT,
										entry.ep.predicate != null ? entry.ep.predicate.shortName : "=");

								if (entry.ep.predicate != null)
								{
									if (entry.ep.predicate.isInterval())
										dtExpValues.setValue("VALUE2", values[1] = entry.ep.getConvertedValue(entry.ep.secondValue, modelMapping.unit));

									// Interchange the "less" and "more" predicates, if the unit was converted and ordering changed
									if (entry.ep.predicate.isOrdering())
										if (!modelMapping.unit.equals(entry.ep.unit) && modelMapping.unit.isOppositeTo(entry.ep.unit))
											dtExpValues.getCurrentRow().addAttachment(QSPRConstants.PREDICATE_ATTACHMENT, entry.ep.predicate.getOpposite());
								}
							}
							else
							{
								dtExpValues.setValue(QSPRConstants.VALUE, values[0] = values[1] = entry.ep.getConvertedAverageValue(modelMapping.unit));
								dtExpValues.getCurrentRow().addAttachment(QSPRConstants.PREDICATE_ATTACHMENT, "="); // We do not handle ranges, always use the "=" predicate
							}

							if(values[0] == null || !Double.isFinite(values[0])) {
								dtExpValues.getCurrentRow().setError("Experimental value: (" + entry.ep.value + ") failed to convert to numerical (" + values[0] +
										") using the specified system of units.");
								continue;
							}

							if(values[1] == null || !Double.isFinite(values[1])) { 
								dtExpValues.getCurrentRow().setError("Experimental value: (" + entry.ep.value+ " - " + entry.ep.secondValue + 
										") failed to convert to numerical (" + values[0] + " - " + values[1] +
										") using the specified system of units.");
								continue;
							}

						}
					}

				if (mapper.getRowsSize() > 1)
					dtExpValues.setValue(QSPRConstants.CLASS, mapper.getClass(entry.ep));
			}catch(Throwable e) {
				dtExpValues.getCurrentRow().setError(e.getMessage());
			}
		}

		return dtExpValues;
	}

	protected Task createTeacherTask() throws Exception
	{
		Task task = new Task(getTaskType(), null, null);
		return task;
	}

	public Task createApplierTask() throws Exception
	{
		Task task = model.template.applierTaskTemplate.getObject();
		return task;
	}

	public void setStatus(String status)
	{
		super.setStatus(status);
		model.detailedStatus = status;
	}

	public void exitAfterPosting()
	{
		exitAfterPosting = true;
		if (model.taskId != null && model.taskId != 0)
			// Task already posted. Exit
			interrupt();
	}

	class MyCalculationClient extends CalculationClient
	{
		public String originalStatus = null;

		public MyCalculationClient()
		{
			super();
			sid = Globals.getClientSID();
			user = Globals.getUsername();
		}

		public void setStatus(String taskStatus)
		{
			super.setStatus(taskStatus);
			String newStatus = originalStatus;
			if (taskStatus != null)
				newStatus += " - " + taskStatus;
			if (!newStatus.equals(ModelProcessor.this.getStatus()))
				ModelProcessor.this.setStatus(newStatus);
		}
	}

	public static ModelProcessor getTeacher(Long modelId)
	{
		for (WrapperThread thread : runningThreads)
			if (thread instanceof ModelProcessor)
			{
				if (modelId.equals(((ModelProcessor) thread).model.id))
					return (ModelProcessor) thread;
			}

		return null;
	}


	private static Logger logger = LogManager.getLogger(ModelProcessor.class);

}
