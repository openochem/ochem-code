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

package qspr.services;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.business.PendingTaskPeer;
import qspr.dao.Repository;
import qspr.entities.BasketEntry;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.Molecule;
import qspr.entities.PendingTask;
import qspr.entities.PendingTask.TaskType;
import qspr.entities.Property;
import qspr.entities.Session;
import qspr.entities.Unit;
import qspr.entities.User;
import qspr.exceptions.CalculationException;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.configurations.ADConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.metaserver.transport.DataReference;
import qspr.metaserver.transport.DataReferenceFactory;
import qspr.modelling.ApplicabilityDomain;
import qspr.modelling.ModelFactory;
import qspr.modelling.ModelProcessor;
import qspr.modelling.ModelStatistics;
import qspr.modelling.PointSelector;
import qspr.modelling.PointStatistics;
import qspr.modelling.SetStatistics;
import qspr.modelling.applier.FeatureNetApplierTaskProcessor;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.applier.ModelApplierTaskProcessor;
import qspr.modelling.applier.Prediction;
import qspr.modelling.applier.PropertyPrediction;
import qspr.toxicity.DispatcherServletWrapper;
import qspr.util.unitconversion.UnitConversion;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.NumericalValueStandardizer;
import com.eadmet.utils.MAILERConstants;
import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;

/**
 * @author vijyant/midnighter/itetko
 * 
 *  A basic web-service for applying a model to get a prediction WSDL available at http://ochem.eu/services/ModelService?wsdl
 */
public class ModelService extends LoginService{
	private static transient final Logger logger = LogManager.getLogger(ModelService.class);
	public boolean useApplierCache = true;

	public String login(String username, String password) {
		return super.login(username, password);
	}

	/**
	 * Required for calculation of consensus models, feature nets, etc.
	 * @param request
	 * @return
	 */
	
	public synchronized ModelResponse postDataReferenceRequest(final DataReferenceRequest request) {
		log("Starting postDataReferenceRequest");
		return new WebServiceWrapper(request.getSessionGUID(), false) 
		{
			@Override
			public void run() throws Exception {
				ThreadScope.get().userSession.disableQuota = true;
				Model model = Repository.model.getByPublicId(request.getModelId());
				model = Globals.initializeAndUnproxy(model);
				if (model == null)
					throw new UserFriendlyException("Unknown model ID: " + request.getModelId());

				log("Use model " + request.getModelId() + " datasize: " + request.getDatasize());

				DataReference reference = DataReferenceFactory.createReference(request.getDataReference(), QSPRConstants.WEBSERVICE_DATABASE);

				FeatureNetApplierTaskProcessor processor = new FeatureNetApplierTaskProcessor(reference, model, 
						PredictionScenario.valueOf(request.getPredictionScenario()), request.getDatasize());
				processor.exitAfterPosting = true;

				processor.wrapped(); // start posting the task
				response.setTaskId(processor.pTask.id);
				response.setMetaserverTaskId(processor.pTask.taskId);
				response.setStatus(QSPRConstants.SUCCESS);
				response.setDetailedStatus(processor.pTask.getDetailedStatus());

				log("Task is posted");

			}
		}.execute();
	}

	/**
	 * Apply model for one sdf file and wait for results
	 * @param modelId
	 * @param sdf
	 * @return
	 */

	public ModelResponse applyModelSingleSDF(final String sessionGUID, final long modelId, final String sdf) {
		log("Starting applyModelSingleSDF");
		String[] sdfs = sdf.split("\\$\\$\\$\\$\n");
		ModelResponse resp = postModelWithSession(null, modelId, sdfs);

		final long taskId = resp.getTaskId();
		if (QSPRConstants.SUCCESS.equals(resp.getStatus())) {
			return new WebServiceWrapper() {
				public void run() throws Exception {
					session = Session.getByGUID(sessionGUID);
					ThreadScope.get().userSession = session;
					
					PendingTask pTask = (PendingTask) Globals.session().get(PendingTask.class, taskId);
					ModelApplier applier = new ModelApplier(pTask);
					Globals.commitAllTransactions();
					while (!applier.isReady()) {
						Thread.sleep(1000);
						Globals.startAllTransactions();
						applier.update();
						Globals.commitAllTransactions();
					}
					Globals.startAllTransactions();
					pTask = (PendingTask) Globals.session().get(PendingTask.class, taskId);
					applier = new ModelApplier(pTask);
					processResults(applier, response, null);
				}
			}.execute();
		} else
			return resp;
	}

	/**
	 * Accept an array of SDFs with a session
	 * @param sessionGUID
	 * @param modelId
	 * @param sdfs
	 * @return
	 */

	public ModelResponse postModelWithSession(final String sessionGUID, final long modelId, final String[] sdfs) 
	{
		log("Starting postModelWithSession");
		PredictionRequest request = new PredictionRequest();
		request.modelId = modelId;
		request.sessionGUID = sessionGUID;
		request.sdfs = sdfs;
		return postModelRequest(request);
	}

	/**
	 * Post a model request using PredictionRequest (session, molecules, sdfs)
	 */

	public synchronized ModelResponse postModelRequest(final PredictionRequest request) {
		log("Starting postModelRequest");
		return new WebServiceWrapper(request.sessionGUID, false) 
		{
			@Override
			public void run() throws Exception {
				if(sessionGUID == null || !sessionGUID.equals(QSPRConstants.ANONYMOUS)) {
					Session sess = Session.getByGUID(sessionGUID);
					if(sess ==  null)throw new UserFriendlyException("Session does not exist: " + sessionGUID);
				}

				ThreadScope.get().userSession.disableQuota = true;
				Model model = Repository.model.getByPublicId(request.getModelId());
				if (model == null)
					throw new UserFriendlyException("Unknown model ID: " + request.modelId);

				ModelApplier applier = new ModelApplier();
				if (request.predictionScenario != null)
					applier.scenario = PredictionScenario.valueOf(request.predictionScenario);
				applier.useCache = useApplierCache;
				applier.addModel(model);

				DataTable tab = new DataTable();
				tab.compressStrings = false;
				tab.addColumn(QSPRConstants.SDF_COLUMN);
				for (String sdf : request.sdfs)
					tab.addRow().setValue(0,sdf);

				applier.defaultTaskPriority = TaskPriority.HIGH; 

				applier.start(tab);

				applier.awaitTasksPosted();
				response.setTaskId(applier.modelTasks.get(0).pTask.id);
				response.setMetaserverTaskId(applier.modelTasks.get(0).pTask.taskId);
				response.setStatus(QSPRConstants.SUCCESS);
				response.setDetailedStatus(applier.modelTasks.get(0).pTask.getDetailedStatus());
			}
		}.execute();
	}

	/*  Deletes a model
	 * @param sessionGUID
	 * @param publicModelID
	 * @return
	 */

	public ModelResponse deleteTask(String sessionGUID, final long taskId) {
		log("Starting deleteTask " + taskId);
		return new WebServiceWrapper(sessionGUID, true) {
			@Override
			public void run() throws Exception {
				PendingTask pTask = (PendingTask) Globals.session().get(PendingTask.class, taskId);
				if (pTask == null)
					throw new UserFriendlyException(Command.UNKNOWN_TASK + taskId);

				checkUser(session.user, pTask.session.user);

				if (pTask.isActive())
					PendingTaskPeer.terminateTaskAsync(pTask.taskId);
				Globals.session().delete(pTask);
			}
		}.execute();
	}

	void checkUser(User session, User owner){
		if(session != null && owner != null && !session.equals(owner))
			throw new UserFriendlyException("Unsufficient privileges to delete task by user " + session + " != " + owner);

		if(session == null && owner != null)
			throw new UserFriendlyException("Unsufficient privileges to delete task by guest user");
	}
	
	
	public ModelResponse getPredictions(final String sessionGUID, final long taskId, final String[] measurementUnit) {
		log("Starting getPredictions for task: " + taskId);
		return new WebServiceWrapper() {
			@Override
			public void run() throws Exception {
				session = Session.getByGUID(sessionGUID);
				ThreadScope.get().userSession = session;

				PendingTask pTask = (PendingTask) Globals.session().get(PendingTask.class, taskId);
				if (pTask == null)
					throw new UserFriendlyException(Command.UNKNOWN_TASK + taskId);

				ModelApplier applier = new ModelApplier(pTask);
				if (applier.isReady())
				{
					if (applier.modelTasks.get(0).pTask.taskId == 0)
						throw new CalculationException("Failed to start model #" + applier.modelTasks.get(0).model.publicId);
					processResults(applier, response, measurementUnit);
				}
				else
				{
					applier.update();
					response.setStatus("pending");
					response.setDetailedStatus(applier.modelTasks.get(0).pTask.getDetailedStatus());
				}
			}
		}.execute();
	}


	public ModelResponse fetchModel(String sessionGUID, final long taskId) {
		log("Starting fetchModel " + taskId);
		return getPredictions(sessionGUID, taskId,null);
	}

	/**
	 * Returns summary of predictions for all molecules in the dataset
	 * 
	 * @param sessionGUID
	 * @param publicModelId
	 * @return
	 */
	public ModelSummary[] getModelSummary(String sessionGUID, Long publicModelId) {
		log("Starting getModelSummary");
		return getModelSummary(sessionGUID, publicModelId, null);
	}


	/**
	 * Calculate model summary for the points matching the provided selector.
	 * If selector is null, all the points are considered
	 * @return
	 */
	private ModelSummary[] getModelSummary(String sessionGUID, Long publicModelId, PointSelector pointSelector)
	{
		log("Starting getModelSummary");
		Globals.startAllTransactions();
		try {
			ThreadScope.get().threadClassLoader = DispatcherServletWrapper.class.getClassLoader();
			Session session = Session.getByGUID(sessionGUID);
			ThreadScope.get().userSession = session;
			if (session == null)
				throw new UserFriendlyException(QSPRConstants.NOSESSION + " " + sessionGUID);

			log("Getting model statistics for model " + publicModelId);
			Model model = Repository.model.getByPublicId(publicModelId);

			int number = 0;

			for (ModelMapping mm : model.modelMappings) {
				ModelStatistics ms = (ModelStatistics) mm.statisticsRecalculated.getObject();
				number += ms.sets.size();
			}

			ModelSummary[] summaries = new ModelSummary[number];
			number = 0;

			/*
			 * For multilearning we can have several models here
			 */

			for (ModelMapping mm : model.modelMappings) {
				ModelStatistics ms = (ModelStatistics) mm.statisticsRecalculated.getObject();
				// TODO: Is this a really correct fix? For the empty statistics values? (NoS)
				if (!ms.containsPredictions) {
					List<PendingTask> pTasks = PendingTask.getByModel(model, TaskType.MODEL_TRAINING);
					if (pTasks.isEmpty())
						throw new UserFriendlyException("Invalid model " + model.publicId);

					log("Getting the task from metaserver");
					Task task = new CalculationClient().getTask(pTasks.get(0).taskId);
					if (!task.isReady())
						throw new UserFriendlyException("The model is still running");
					if (task.isError())
						throw new UserFriendlyException("The model has failed: " + task.getDetailedStatus());
					ModelProcessor processor = ModelFactory.getProcessor(model.template);
					processor.model = model;
					processor.onTaskReceived(task);
					processor.saveModel();

					// Reload the statistics
					ms = (ModelStatistics) mm.statisticsRecalculated.getObject();
				}

				ms.recalculateStatistics(mm, pointSelector);
				int set = 0;
				for (SetStatistics ss : ms.sets) {

					summaries[number] = new ModelSummary();

					summaries[number].setModelName(model.name);

					summaries[number].setPropertyName(mm.property.getName());
					summaries[number].setValidationSetName(set == 0 ? QSPRConstants.TRAINING : QSPRConstants.VALIDATION + set);

					summaries[number].setRmse(NumericalValueStandardizer.getSignificantDigitsDouble(ss.rmse, NumericalValueStandardizer.SIGNIFICANT_DIGITS));
					summaries[number].setRmseConfidence(NumericalValueStandardizer.getSignificantDigitsDouble(ss.rmseStd,NumericalValueStandardizer.SIGNIFICANT_DIGITS));
					summaries[number].setMae(NumericalValueStandardizer.getSignificantDigitsDouble(ss.mae,NumericalValueStandardizer.SIGNIFICANT_DIGITS));
					summaries[number].setMaeConfidence(NumericalValueStandardizer.getSignificantDigitsDouble(ss.maeStd,NumericalValueStandardizer.SIGNIFICANT_DIGITS));
					summaries[number].setR2(NumericalValueStandardizer.getSignificantDigitsDouble(ss.r2,NumericalValueStandardizer.SIGNIFICANT_DIGITS));
					summaries[number].setR2Confidence(NumericalValueStandardizer.getSignificantDigitsDouble(ss.r2,NumericalValueStandardizer.SIGNIFICANT_DIGITS));
					summaries[number].setQ2(NumericalValueStandardizer.getSignificantDigitsDouble(ss.q2,NumericalValueStandardizer.SIGNIFICANT_DIGITS));
					summaries[number].setQ2Confidence(NumericalValueStandardizer.getSignificantDigitsDouble(ss.q2,NumericalValueStandardizer.SIGNIFICANT_DIGITS));
					summaries[number].setN(ss.n);

					if (ss.classificationSummary != null) {
						summaries[number].setAccuracy(NumericalValueStandardizer.getSignificantDigitsDouble(ss.classificationSummary.accuracyTotal.getValue(),NumericalValueStandardizer.SIGNIFICANT_DIGITS));
						summaries[number].setAccuracyConfidence(NumericalValueStandardizer.getSignificantDigitsDouble(ss.classificationSummary.accuracyTotal.getStd(),NumericalValueStandardizer.SIGNIFICANT_DIGITS));
						summaries[number].setAccuracyBalanced(NumericalValueStandardizer.getSignificantDigitsDouble(ss.classificationSummary.accuracyBalanced.getValue(),NumericalValueStandardizer.SIGNIFICANT_DIGITS));
						summaries[number].setAccuracyBalancedConfidence(NumericalValueStandardizer.getSignificantDigitsDouble(ss.classificationSummary.accuracyBalanced.getStd(),NumericalValueStandardizer.SIGNIFICANT_DIGITS));
						summaries[number].setMcc(NumericalValueStandardizer.getSignificantDigitsDouble(ss.classificationSummary.mcc.getValue(),NumericalValueStandardizer.SIGNIFICANT_DIGITS));
						summaries[number].setMccConfidence(NumericalValueStandardizer.getSignificantDigitsDouble(ss.classificationSummary.mcc.getStd(),NumericalValueStandardizer.SIGNIFICANT_DIGITS));
						summaries[number].setTP(ss.classificationSummary.tp);
						summaries[number].setFP(ss.classificationSummary.fp);
						summaries[number].setTN(ss.classificationSummary.tn);
						summaries[number].setFN(ss.classificationSummary.fn);
					}
					number++;
					set++;
				}
			}

			Globals.commitAllTransactions();

			return summaries;

		} catch (RuntimeException e) {
			Globals.rollbackAllTransactions();
			throw e;
		} catch (Exception e) {
			Globals.rollbackAllTransactions();
			e.printStackTrace();
			throw new UserFriendlyException(e.getMessage());
		} finally {
			ThreadScope.reset();
		}
	}

	/**
	 * Provides predictions
	 * @param sessionGUID
	 * @param publicModelId
	 * @return
	 */

	public PropertyPrediction[] getModelPredictions(String sessionGUID, Long publicModelId) {
		log("Starting getModelPredictions");
		Session session = null;

		Globals.startAllTransactions();
		ThreadScope.get().threadClassLoader = DispatcherServletWrapper.class.getClassLoader();
		try {
			session = Session.getByGUID(sessionGUID);
			ThreadScope.get().userSession = session;
			if (session == null)
				throw new UserFriendlyException(QSPRConstants.NOSESSION + " " + sessionGUID);

			log("Getting model statistics for model " + publicModelId);
			Model model = Repository.model.getByPublicId(publicModelId);

			List<PropertyPrediction> ppList = new ArrayList<PropertyPrediction>();

			HashMap<Long, Integer> epToMol = new HashMap<Long, Integer>();
			for (ModelMapping mm : model.modelMappings) {
				ModelStatistics ms = (ModelStatistics) mm.statisticsRecalculated.getObject();
				for (SetStatistics ss : ms.sets)
					for (PointStatistics ps : ss.points)
						epToMol.put(ps.id, 0);
			}

			log("Getting molecule IDs");
			@SuppressWarnings("unchecked")
			List<Object[]> rows = Globals.session().createCriteria(ExperimentalProperty.class).createAlias("molecule", "mol")
			.createAlias("mol.mapping2", "mp2").add(Restrictions.in("id", epToMol.keySet()))
			.setProjection(Projections.projectionList().add(Projections.groupProperty("mp2.id")).add(Projections.groupProperty("id"))).list();

			for (Object[] row : rows) {
				epToMol.put((Long) row[1], (Integer) row[0]);
			}

			log("Preparing statistics");
			for (ModelMapping mm : model.modelMappings) {
				ModelStatistics ms = (ModelStatistics) mm.statisticsRecalculated.getObject();
				if (!ms.containsPredictions) {
					List<PendingTask> pTasks = PendingTask.getByModel(model, TaskType.MODEL_TRAINING);
					if (pTasks.isEmpty())
						throw new UserFriendlyException("Invalid model " + model.publicId);

					log("Getting the task from metaserver");
					Task task = new CalculationClient().getTask(pTasks.get(0).taskId);
					if (!task.isReady())
						throw new UserFriendlyException("The model is still running");
					if (task.isError())
						throw new UserFriendlyException("The model has failed: " + task.getDetailedStatus());
					ModelProcessor processor = ModelFactory.getProcessor(model.template);
					processor.model = model;
					processor.onTaskReceived(task);
					processor.saveModel();

					// Reload the statistics
					ms = (ModelStatistics) mm.statisticsRecalculated.getObject();
				}

				ApplicabilityDomain ad = null;
				if (ms.sets.get(0).distancesToModel.size() > 0) {

					int dmNum = 0;
					String dmName = ms.sets.get(0).distancesToModel.get(dmNum);

					log("Calculating applicability domain: " + dmName);

					ad = new ApplicabilityDomain();
					ADConfiguration adConfiguration = new ADConfiguration();
					ad.setModel(mm, dmName, adConfiguration);
				}

				for (SetStatistics ss : ms.sets) {
					for (PointStatistics ps : ss.points) {
						PropertyPrediction pp = new PropertyPrediction();
						if (ps.error == null) {
							pp.setValue(ps.predicted); // FIXME: Multiple classification case is not supported so far.  IVT: changed to raw value by request of Ahmed
							pp.setRealValue(ps.real);

							if (mm.property.isQualitative()) {
								pp.setRealValueString(model.attachment.getObject().getOptionFromPrediction(ps.real, mm.property).name);
								pp.setPredictedValueString(model.attachment.getObject().getOptionFromPrediction(ps.predicted, mm.property).name);
							}

							pp.setProperty(mm.property.getName());
							pp.setMoleculeId("M" + epToMol.get(ps.id));
							if (ad != null)
								pp.setAccuracy(ad.getPredictedError(ps.distancesToModel.get(0)));

							if (ad != null && ps.distancesToModel != null)
								pp.setDm(ps.distancesToModel.get(0));
						} else
							pp.setError(ps.error);

						ppList.add(pp);
					}

				}
			}

			Globals.commitAllTransactions();

			PropertyPrediction[] result = new PropertyPrediction[ppList.size()];
			for (int i = 0; i < ppList.size(); i++) {
				result[i] = ppList.get(i);
			}

			return result;

		} catch (RuntimeException e) {
			Globals.rollbackAllTransactions();
			throw e;
		} catch (IOException e) {
			Globals.rollbackAllTransactions();
			e.printStackTrace();
			throw new UserFriendlyException("Metaserver is not available");
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			throw new UserFriendlyException(e.getMessage());
		} catch (Exception e) {
			Globals.rollbackAllTransactions();
			e.printStackTrace();
			throw new UserFriendlyException(e.getMessage());
		} finally {
			ThreadScope.reset();
		}
	}

	private Unit[] convertUnits(String[] measurementUnit,  Model model){

		if(measurementUnit == null || measurementUnit.length == 0) return null;

		// Retrieve the units by their names
		Unit[] units = null;
		if (measurementUnit != null && measurementUnit.length > 0)
		{
			units = new Unit[measurementUnit.length];
			for (int i = 0; i < units.length; i++) 
				if (measurementUnit[i] != null && !"".equals(measurementUnit[i]))
					units[i] = Unit.getByNameAndCategory(measurementUnit[i], model.modelMappings.get(i).property.unitCategory.name, false);

			for(int pNum = 0; pNum < measurementUnit.length; pNum++) {
				ModelMapping mm = model.modelMappings.get(pNum);
				if(units[pNum]  != null && mm.unit.id != units[pNum].id)return units;
			}
		}

		return null;
	}

	private void processResults(ModelApplier applier, ModelResponse response, String[] measurementUnit) throws Exception {
		ModelApplierTaskProcessor mTask = applier.modelTasks.get(0);

		if (mTask.isError())
			throw new UserFriendlyException("Calculation failed "+mTask.getStatus());

		Model model = mTask.model = (Model) Globals.session().get(Model.class, mTask.model.id);  // refetch model
		DataTable dtResult = mTask.wndResult.ports.get(0);

		Unit[] units = convertUnits(measurementUnit, model);

		// Initialize applicability domains
		mTask.resetApplicabilityDomains(null);
		List<BasketEntry> entries = applier.entries();

		if(entries.size() != dtResult.getRowsSize()){  // attachment was not saved for this task -- we retrieve the parent task which have it for consensus model
			String msg = "For " + mTask.pTask.taskId+ " results = "+dtResult.getRowsSize()+" in attachment " + entries.size();
			Mailer.postMailSafely(new Email(MAILERConstants.EMAIL_ADMIN, "Exception in ModelService for task " + mTask.pTask.taskId, msg));
			throw new Exception("Unexpected error with calculation: contact us at "+MAILERConstants.EMAIL_OCHEM+ " and provide details of calculations. Thank you!");
		}

		Prediction[] predictions = new Prediction[dtResult.getRowsSize()];
		response.setPredictions(predictions);

		boolean hasGood = false;
		dtResult.reset();
		while (dtResult.nextRow()) {
			Prediction prediction = predictions[dtResult.currentRow]  = new Prediction();

			Molecule molecule = entries.get(dtResult.currentRow).ep.molecule;
			// Set the molecule identity info

			prediction.setSmiles((String)dtResult.getRow(dtResult.currentRow).attachments.get(QSPRConstants.SMILES_ATTACHMENT)); 
			String inchi = (String)dtResult.getRow(dtResult.currentRow).attachments.get(QSPRConstants.INCHIKEYS);
			if(inchi != null && inchi.contains("-"))
				prediction.setInChIKey(inchi); 

			if(molecule != null && !dtResult.getCurrentRow().isError()) {
				hasGood = true;
				//prediction.setDepictionID(molecule.id);
				prediction.setMoleculeID(molecule.mapping2.id);
				response.setStatus(QSPRConstants.SUCCESS);
				int pNum = 0;
				List<PropertyPrediction> propPredictions = new ArrayList<PropertyPrediction>();
				for (String c : dtResult.getColumns()) {

					if (c.toLowerCase().startsWith("result")) {
						ModelMapping mm = model.modelMappings.get(pNum);
						Property property = mm.property;
						ApplicabilityDomain ad = mTask.getApplicabilityDomain(null, pNum);
						PropertyPrediction pPred = new PropertyPrediction();
						double predictionValue = ((Double) dtResult.getValue(c)).doubleValue(); 
						pPred.setValue(NumericalValueStandardizer.getSignificantDigitsDouble(predictionValue, NumericalValueStandardizer.SIGNIFICANT_DIGITS));
						pPred.setPredictedValueString(NumericalValueStandardizer.getSignificantDigitsStr(predictionValue, NumericalValueStandardizer.SIGNIFICANT_DIGITS));
						pPred.setUnit(model.modelMappings.get(pNum).unit.getName());
						pPred.setProperty(property.getName());
						if (property.isQualitative())
							pPred.setPredictedValueString(model.attachment.getObject().getOptionFromPrediction(pPred.getValue(), property).name);
						else
						{
							// Convert the units, if required
							if (units != null && units.length > pNum && units[pNum] != null)
								pPred.setValue(UnitConversion.convert(pPred.getValue(), mm.unit, units[pNum], molecule.molWeight));
						}

						if (ad != null) {
							Double error = ad.getPredictedError(dtResult.currentRow);
							pPred.setAccuracy(error != null ? NumericalValueStandardizer.getSignificantDigitsDouble(error, NumericalValueStandardizer.SIGNIFICANT_DIGITS) : 0.0d);
							Double dm = ad.getPredictedDM(dtResult.currentRow);
							pPred.setDm(dm != null ? NumericalValueStandardizer.getSignificantDigitsDouble(dm, NumericalValueStandardizer.SIGNIFICANT_DIGITS) : 0.0d);
							if (dm != null)
								pPred.setInsideAD(dm <= ad.dmThreshold);
						}
						propPredictions.add(pPred);
						pPred.addRealValues(property, molecule, mm.unit,model.trainingSet);

						pNum++;
					}
				}

				PropertyPrediction[] predArray = new PropertyPrediction[propPredictions.size()];
				for (int i = 0; i < propPredictions.size(); i++)
					predArray[i] = propPredictions.get(i);
				prediction.setPredictions(predArray);
			}else{
				String error = dtResult.getCurrentRow().detailedStatus;
				if(error == null ||  error.length() == 0) error = QSPRConstants.ERROR_STATUS;
				prediction.setError(error.trim());
				log(prediction.getError());
			} 
		}
		if(!hasGood)response.setStatus("failed");
		else
			response.setStatus(QSPRConstants.SUCCESS);
		response.setModelDescriptionUrl(OCHEMConfiguration.rootHost+"/model/" + applier.modelTasks.get(0).model.publicId);
	}

	private void log(String msg) {
		if(msg !=null && Globals.userSession() != null && Globals.userSession().user != null && Globals.userSession().user.login != null) msg += " by user " + Globals.userSession().user.login;
		System.out.println("message "+msg);
		logger.info("[WebService] " + msg);
	}
}

