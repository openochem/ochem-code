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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.dao.Repository;
import qspr.entities.Attachment;
import qspr.entities.AttachmentSource;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.Molecule;
import qspr.entities.PendingTask;
import qspr.entities.Attachment.AttachmentType;
import qspr.export.ExportableMolecule;
import qspr.export.ExportableSetConfiguration;
import qspr.interfaces.ProvidedConditions;
import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.util.ShortCondition;
import qspr.workflow.utils.QSPRConstants;
import qspr.modelling.ApplicabilityDomain;
import qspr.modelling.CDSModelProcessor;
import qspr.modelling.CompoundsProvider;
import qspr.modelling.ModelApplierAttachment;
import qspr.modelling.ModelApplierCache;
import qspr.modelling.ModelProcessor;
import qspr.modelling.configurations.ExternalCondition;
import qspr.util.ExportThread;
import qspr.util.MoleculePeer;
import qspr.util.WrapperThread;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.NumericalValueStandardizer;

public class ModelApplier
{
	private static transient final Logger logger = LogManager.getLogger(ModelApplier.class);

	public int BATCH_SIZE_STORED_MOLECULE = QSPRConstants.WEBSERVICE_SIZE_STORED_MOLECULE;

	
	public PredictionScenario scenario = PredictionScenario.PREDICTION_ONLY;
	public ModelProcessor teacher = new CDSModelProcessor();
	public ArrayList<ModelApplierTaskProcessor> modelTasks = new ArrayList<ModelApplierTaskProcessor>();
	public int repostSize = QSPRConstants.MODEL_REPOST_SIZE; 
	// data can be supplied as Basket
	public CompoundsProvider compoundsProvider = new CompoundsProvider();
	public Integer defaultTaskPriority;
	public ConditionSet defaultConditions;

	public List<Integer> order;

	private String resultOrdering;
	private boolean ascending;

	public boolean refreshModelsFromDB = true;
	public boolean useCache = true;
	public boolean forceUpdateDescriptorCache = false;

	public int taskDebugLevel = 0;

	private ModelApplierAttachment attachment;

	public Integer indexOf(Model m)
	{
		for (int i=0; i < modelTasks.size(); i++) 
		{
			ModelApplierTaskProcessor elem = modelTasks.get(i);
			if (elem.model.id.equals(m.id))
				return Integer.valueOf(i);
		}
		return null;
	}

	public void addModel(Model m)
	{
		ModelApplierTaskProcessor mt = new ModelApplierTaskProcessor(this, m);
		modelTasks.add(mt);
	}

	public void addModelApplierTask(ModelApplierTaskProcessor mt) {
		modelTasks.add(mt);
	}

	public void removeModel(Model m)
	{
		Integer i = indexOf(m);
		if (i != null)
			modelTasks.remove(i.intValue());
	}

	public void start(DataTable molecules) throws Exception
	{
		Long val; 
		molecules.reset();

		for (int i=1; molecules.nextRow(); i++){

			if (++i % 1000 == 0)
				Globals.restartAllTransactions(true); // each 1000 rows

			// we have record!
			if((val = (Long)molecules.getCurrentRow().getAttachment(QSPRConstants.RECORD_ID_ATTACHMENT)) !=null){
				ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, val);
				entries().add(new BasketEntry(ep));
				continue;

			}

			// we have Molecules
			ExperimentalProperty ep = new ExperimentalProperty();

			if((val = (Long)molecules.getCurrentRow().getAttachment(QSPRConstants.MOLECULE_ID_OCHEM)) !=null){
				Molecule mol = Repository.molecule.getMolecule(val);
				ep.molecule = mol;
				entries().add(new BasketEntry(ep));
				continue;
			}

			// we have only SDF
			String sdf = (String)molecules.getValue();

			if (sdf == null)
				ep.molecule = Repository.molecule.getEmptyMolecule();
			else if (sdf.matches("M[0-9]+"))
				ep.molecule = MoleculePeer.getByID2(sdf);
			else
			{
				Molecule mol = null;

				// Heuristics. Fetch the molecules from DB only for small data-sets
				if (molecules.getRowsSize() > BATCH_SIZE_STORED_MOLECULE) {
					// Create a molecule stub so that its not saved in the OCHEM DB
					mol = Molecule.getStub();
					mol.setData(sdf);
					mol.mapping2.id = -i;
				}
				else
				{
					try
					{
						mol = MoleculePeer.getMolecule(sdf);
					}
					catch (Exception e) {
						mol = Molecule.getStub();
						mol.setData("AAA");
						mol.error = e.getMessage();
						mol.mapping2.id = -i;
					}
				}
				ep.molecule = mol;

			}

			entries().add(new BasketEntry(ep));
		}

		Globals.restartAllTransactions(true); // final part

		start();
	}

	/**
	 * Starting a thread with the task
	 */

	public void start()
	{
		logger.info("Starting predictions by " + modelTasks.size() + " models");
		attachment = null;
		for (ModelApplierTaskProcessor modelTask : modelTasks)
		{
			modelTask.taskDescription = getConditionsDescription();
			modelTask.exitAfterPosting = true;
			modelTask.savePendingTaskEarly = true;
			modelTask.start();
		}
	}

	/**
	 * A shorthand for accessing basket entries
	 * @return
	 */
	public List<BasketEntry> entries()
	{
		return compoundsProvider.basket.entries;
	}

	/**
	 * Start and wait till completion. This is a blocking function that can take unlimited time to execute.
	 */
	public void startAndWait() throws Exception
	{
		start();
		awaitTasksFinished();
	}

	public void update() throws Exception
	{
		for (ModelApplierTaskProcessor modelTask : modelTasks)
			modelTask.update();
	}

	public void awaitTasksPosted() throws Exception
	{
		boolean pending = true;

		while (pending)
		{
			pending = false;
			for (ModelApplierTaskProcessor modelTask : modelTasks)
				if (modelTask.pTask == null || modelTask.pTask.id == null || modelTask.isRunning())
					if (modelTask.isError())
						throw modelTask.exception;
					else
						pending = true;
			if (pending)
				Thread.sleep(200);
		}

		logger.info("Task has been posted "+modelTasks.get(0).pTask);
	}

	public void awaitTasksFinished() throws Exception {
		int timeout = 100;
		while (!isReady())
		{
			Thread.sleep(timeout);
			update();
			timeout = 1000;
		}
	}

	/**
	 * Gets non-redundant list of experimental conditions required by the selected models
	 * @return
	 */
	public List<ExternalCondition> getNecessaryConditions()
	{
		List<ExternalCondition> eDescs = new ArrayList<ExternalCondition>();
		Set<Long> conditions = new HashSet<Long>();

		for (ModelApplierTaskProcessor mTask : modelTasks) 
		{
			mTask.initialiseModel();
			if (mTask.model.attachment.getObject().configuration instanceof ProvidedConditions &&
					((ProvidedConditions) mTask.model.attachment.getObject().configuration).hasConditions()) {

				List<ShortCondition> eD = ((ProvidedConditions) mTask.model.attachment.getObject().configuration).getConditions();
				for(ShortCondition e :eD) {
					if(conditions.contains(e.id))continue;
					conditions.add(e.id);
					eDescs.add(new ExternalCondition(e));
				}
			}
		}

		return eDescs;
	}

	public String getPrintableStatus()
	{
		int ready = 0, error = 0, running = 0;
		for (ModelApplierTaskProcessor modelTask : modelTasks)
		{
			if (modelTask.isReady())
				if (modelTask.wndResult != null)
					ready++;
				else
					error++;
			else
				running++;
		}
		return ""+ready+" ready, "+running+" running, "+error+" failed";
	}

	public boolean isReady()
	{
		for (ModelApplierTaskProcessor mt : modelTasks)
			if (!mt.isReady())
				return false;
		return true;
	}

	public boolean isError()
	{
		for (ModelApplierTaskProcessor mt : modelTasks)
			if (mt.isError())
				return true;
		return false;
	}

	public String getErrorMessage() {
		List<String> errors = new ArrayList<String>();
		for (ModelApplierTaskProcessor mt : modelTasks)
			if (mt.isError())
				errors.add(mt.getStatus());
		return StringUtils.join(errors, "\n");
	}

	public void setResultOrdering(String resultOrdering, boolean ascending) throws Exception
	{
		if (resultOrdering == null)
			resultOrdering = "";


		final int sgn = ascending ? 1 : -1;

		if (resultOrdering.equals(this.resultOrdering) && this.ascending == ascending)
			return;

		if ("".equals(resultOrdering))
			order = null;
		else
		{
			order = new ArrayList<Integer>();
			for (int i = 0; i < compoundsProvider.basket.entries.size(); i++)
				order.add(i);

			if ("dm".equals(resultOrdering))
			{
				ModelApplierTaskProcessor mt = null;
				for (int i = 0; i < modelTasks.size(); i++)
					if ((mt = modelTasks.get(i)).getApplicabilityDomain(null, 0) != null)
						break;
				if (mt.getApplicabilityDomain(null, 0) != null)
				{
					final WorkflowNodeData wndResult = mt.wndResult;
					final ApplicabilityDomain ad = mt.getApplicabilityDomain(null, 0);
					Collections.sort(order, new Comparator<Integer>(){
						public int compare(Integer o1, Integer o2) 
						{
							Double dm1 = ad.getPredictedDM(o1);
							Double dm2 = ad.getPredictedDM(o2);

							if (wndResult.ports.get(0).getRow(o1).isError())
								return 1;

							if (wndResult.ports.get(0).getRow(o2).isError())
								return -1;

							if (dm1 == null)
								return -1;
							if (dm2 == null)
								return 1;

							if (dm1.equals(dm2))
								return 0;
							else
								return -sgn * (dm1 > dm2 ? 1 : -1);
						}
					});	
				}
			}
			else if ("prediction".equals(resultOrdering))
			{
				final ModelApplierTaskProcessor mt = modelTasks.get(0);
				Collections.sort(order, new Comparator<Integer>(){
					public int compare(Integer o1, Integer o2) 
					{
						Double val1 = (Double) mt.wndResult.ports.get(0).getValue(o1, 0);
						Double val2 = (Double) mt.wndResult.ports.get(0).getValue(o2, 0);
						if (val1 == null)
							return -1;
						if (val2 == null)
							return 1;
						if (val1.equals(val2))
							return 0;
						else
							return sgn * (val1 > val2 ? 1 : -1);
					}
				});	
			}

		}

		this.resultOrdering = resultOrdering;
		this.ascending = ascending;
	}

	public void reset()
	{
		compoundsProvider.basket = new Basket();
		for (ModelApplierTaskProcessor mt  : modelTasks) 
		{
			mt.wndResult = null;
		}
	}

	protected synchronized ModelApplierAttachment getAttachment()
	{
		if (attachment != null)
			return attachment;
		attachment = new ModelApplierAttachment();
		attachment.setWorkData(compoundsProvider.getBasket());

		return attachment;
	}

	public ModelApplier()
	{

	}

	/**
	 * Restore ModelAIpplier from attachments
	 * @param pt
	 * @throws IOException
	 * @throws ClassNotFoundException
	 * @throws Exception
	 */

	public ModelApplier(PendingTask pt) throws IOException, ClassNotFoundException, Exception
	{
		addModel(pt.model);
		//modelTasks.get(0).taskId = pt.taskId;
		ModelApplierTaskProcessor processor = modelTasks.get(0);

		processor.pTask = pt;
		processor.setDescription = pt.setDescription;

		if (pt.taskId == 0)
			if (WrapperThread.getPendingTaskThread(pt.id) != null)
				modelTasks.set(0, (ModelApplierTaskProcessor) WrapperThread.getPendingTaskThread(pt.id));

		if (pt.attachment != null)
			attachment = (ModelApplierAttachment) pt.attachment.getObject();

		if (attachment != null)
		{
			compoundsProvider.basket = attachment.getWorkData(false);

			if (attachment.cache != null)
			{
				attachment.cache.model = pt.model;
				attachment.cache.basket = compoundsProvider.basket;
				processor.setCache(attachment.cache);
			}
		}

		if (pt.taskId == -1)
		{
			// Its not a real task, everything was cached
			if(pt.readyTask != null){
				Task task = pt.retrieveTask(true);
				processor.onTaskReceived(task);
			}else{
				//To avoid very long access time, we save it as a normal pseudo - Task
				processor.onTaskReceived(null);
				Task task = new Task();
				task.id = -1; task.status = Task.READY;
				task.setResult(processor.wndResult); // updating all results with new ones
				attachment.cache = null;
				pt.attachment.setObject(attachment, AttachmentType.SERIALIZABLE);
				pt.readyTask = new Attachment<Task>(task, AttachmentType.SERIALIZABLE, AttachmentSource.LocalPendingTask);
				Globals.session().saveOrUpdate(pt);

			}
		}
		else if (pt.taskId != 0)
		{
			Task task = pt.retrieveTask(true);
			if (task != null){
				processor.onTaskReceived(task);
				if(pt.readyTask == null){
					attachment.cache = null;
					pt.attachment.setObject(attachment, AttachmentType.SERIALIZABLE);
					pt.readyTask = new Attachment<Task>(task, AttachmentType.SERIALIZABLE, AttachmentSource.LocalPendingTask);
					Globals.session().saveOrUpdate(pt);
				}
			}
		}
	}

	protected void fetchRecordsFromDatabase()
	{
		logger.info("Fetching records from DB...");
		List<Long> idsToFetch = new ArrayList<Long>();

		for (int batchStart = 0, i = 0; i < compoundsProvider.basket.entries.size(); i++)
		{
			BasketEntry be = compoundsProvider.basket.entries.get(i);
			if (be.ep != null && be.ep.property == null)
				idsToFetch.add(be.ep.id);

			// Load records from DB in batches of 500 records. TODO: Consider restarting transaction (needs careful debugging) / Midnighter on Jun 27, 2011
			if (idsToFetch.size() >= CompoundsProvider.NUMBER_OF_RECORDS_TO_FETCH || i == compoundsProvider.basket.entries.size() - 1)
			{
				if (!idsToFetch.isEmpty())
					Globals.session().createCriteria(ExperimentalProperty.class).add(Restrictions.in("id", idsToFetch)).list();
				for (int k = batchStart; k <= i; k++)
				{
					BasketEntry be2 = compoundsProvider.basket.entries.get(k);
					if (be2.ep != null && be2.ep.property == null)
						be2.ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, be2.ep.id);
					batchStart = i + 1;
					idsToFetch.clear();
				}

				logger.info("" + i + " records fetched");
			}
		}
	}

	public ExportThread getExportThread(ExportableSetConfiguration conf, String format)
	{
		ExportThread eThread = new ExportThread(format, conf)
		{
			@SuppressWarnings({ "unchecked"})
			@Override
			public void generateData() throws Exception
			{
				int BATCH_SIZE = 200;
				int baseIndex = 0;

				eData.setDescriptors(modelTasks.get(0).dtDescriptors);

				Map<Long, Model> modelCache = new HashMap<Long, Model>();
				for (int j = 0; j < modelTasks.size(); j++)
				{
					ModelApplierTaskProcessor mt = modelTasks.get(j);
					mt.initialiseModel();
					modelCache.put(mt.model.id, mt.model);
				}

				while (baseIndex < compoundsProvider.basket.entries.size())
				{
					setStatus("Prepared " + eData.exportableMolecules.size() + " predictions out of " + compoundsProvider.basket.entries.size());
					if (baseIndex % 1000 == 0 && baseIndex > 0)
					{
						Globals.restartAllTransactions(true);
						Runtime.getRuntime().gc();
					}

					List<Long> ids = new ArrayList<Long>();

					for (int i = 0; i < BATCH_SIZE; i++)
						if (baseIndex + i < compoundsProvider.basket.entries.size())
							if (compoundsProvider.basket.entries.get(baseIndex + i).ep.id != null)
								ids.add(compoundsProvider.basket.entries.get(baseIndex + i).ep.id);

					List<ExperimentalProperty> eps = null;

					if (ids.size() > 0)
						eps = Globals.session().createCriteria(ExperimentalProperty.class).add(Restrictions.in("id", ids)).list();

					Map<Long, ExperimentalProperty> epCache = new HashMap<Long, ExperimentalProperty>();

					if (eps != null)
						for (int i = 0; i < eps.size(); i++)
							epCache.put(eps.get(i).id, eps.get(i));

					for (int i = 0; i < BATCH_SIZE; i++)
					{
						if (baseIndex + i >= compoundsProvider.basket.entries.size())
							continue;

						ExperimentalProperty ep = compoundsProvider.basket.entries.get(baseIndex + i).ep; // Observed

						ExportableMolecule eMol = new ExportableMolecule();
						eData.addMolecule(eMol);

						// value
						if (ep.id != null)
						{
							ep = epCache.get(ep.id);
							if(ep == null) {
								eMol.error = "This record has been deleted and cannot be exported.";
								continue;
							}
							eMol.setExperimentalProperty(ep);
						}
						else
						{ //assuming that all records without experimental properties are uploaded; that is wrong! 
							//TODO add a better control of download of data from TAGs
							eMol.setMolecule(ep.molecule);
							eMol.setMoleculeNames(ep.moleculenames);
							eMol.ownRecord = true;
							if(ep.other != null) eMol.comments = ep.other;
						}

						if (modelTasks.get(0).dtDescriptors != null)
							eMol.descriptors = modelTasks.get(0).dtDescriptors.getRow(baseIndex + i);

						for (int j = 0; j < modelTasks.size(); j++)
						{
							// int column = colheader + mbData.modelTasks.size() + 1;
							ModelApplierTaskProcessor mt = modelTasks.get(j);
							mt.model = modelCache.get(mt.model.id);
							DataTable dtResult = null;
							Model m = mt.model;
							boolean isMultilearning = m.modelMappings.size() > 1;
							String result = mt.wndResult == null ? "error" : mt.wndResult.ports.get(0).getRow(baseIndex + i).status;
							boolean errorRow = result != null && result.equals("error");
							if (!errorRow)
								dtResult = mt.wndResult.ports.get(0);
							for (int mmNum = 0; mmNum < m.modelMappings.size(); mmNum++)
							{
								ApplicabilityDomain ad = mt.getApplicabilityDomain(null, mmNum);
								ModelMapping mm = m.modelMappings.get(mmNum);
								if (errorRow)
									eMol.error = mt.wndResult == null ? mt.getStatus() : mt.wndResult.ports.get(0).getRow(baseIndex + i).detailedStatus;
								else
								{
									// A valid result, no error
									// 1) Prediction
									Double val;
									if (isMultilearning)
										val = (Double) dtResult.getValue(baseIndex + i, QSPRConstants.PREDICTION_RESULT_COLUMN + mmNum);
									else
										val = (Double) dtResult.getValue(baseIndex + i, 0);
									val = NumericalValueStandardizer.getSignificantDigitsDouble(val, NumericalValueStandardizer.SIGNIFICANT_DIGITS);

									if (mm.property.isQualitative())
										eMol.setPrediction(mm, m.attachment.getObject().getOptionFromPrediction(val, mm.property).name, val);
									else
										eMol.setPrediction(mm, val);

									// 2) Estimated accuracy
									if (mt.getApplicabilityDomain(null, 0) != null)
									{
										Double pe = Math.rint(ad.getPredictedError(baseIndex + i) * 100) / 100;
										eMol.setAccuracy(mm, pe);
										if (ad.dmThreshold != null && ad.getPredictedDM(baseIndex + i) > ad.dmThreshold)
											eMol.ousideOfAD = true;
									}

									// DMs
									String dmPrefix = "DM:";
									if (m.modelMappings.size() > 1)
										dmPrefix = "DM" + mmNum + ":";
									for (int k = 0; k < dtResult.getColumnsSize(); k++)
									{
										if (dtResult.getColumn(k).startsWith(dmPrefix))
										{
											String dmName = dtResult.getColumn(k).substring(dmPrefix.length());
											eMol.setDM(mm, dmName, Math.rint((Double) dtResult.getValue(baseIndex + i, k) * 100) / 100);
										}
									}


								}
							}
						}
					}
					baseIndex += BATCH_SIZE;
					logger.info(MemoryUtils.memorySummary());
				}
				setFileName("ModelPredictions_" + eData.exportableMolecules.size() + "_compounds");
			}

		};
		return eThread;
	}

	public static ModelApplier getFromTask(long pendingTaskId) throws Exception {
		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		if (applier != null && applier.modelTasks.get(0).pTask != null && applier.modelTasks.get(0).pTask.id.equals(pendingTaskId))
			return applier;
		else
			return new ModelApplier(PendingTask.getById(pendingTaskId));
	}

	public static void clearCachedPredictions(Long publicId) {
		System.out.println("Deleting cache for model " + publicId);
		ModelApplierCache.clearCachedPredictions(publicId);
	}

	public String getConditionsDescription() {
		if(defaultConditions == null) return null;
		return "Used conditions: " + defaultConditions;
	}

}