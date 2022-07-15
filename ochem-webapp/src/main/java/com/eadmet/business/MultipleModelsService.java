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

package com.eadmet.business;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.business.PendingTaskPeer;
import qspr.dao.ChemInfEngine;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.ModelConfigurationTemplate;
import qspr.entities.ModelConfigurationTemplate.TemplateType;
import qspr.entities.ModelMapping;
import qspr.entities.PendingTask;
import qspr.entities.PendingTask.TaskType;
import qspr.entities.Property;
import qspr.entities.Unit;
import qspr.export.ExportableModel;
import qspr.frontend.MultipleModelsData;
import qspr.frontend.MultipleModelsData.ModelData;
import qspr.frontend.MultipleModelsData.ModelType;
import qspr.frontend.MultipleModelsData.SetStats;
import qspr.interfaces.ProvidedConditions;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.configurations.BaggingConfiguration;
import qspr.metaserver.configurations.CompressedObject;
import qspr.metaserver.configurations.ConsensusModelConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration.MixturesProcessing;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.NoDescriptors;
import qspr.metaserver.configurations.StandartizationOptions;
import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.metaserver.configurations.ValidationConfiguration;
import qspr.modelling.CrossModelGenerator;
import qspr.modelling.ModelFactory;
import qspr.modelling.ModelProcessor;
import qspr.modelling.ModelStatistics;
import qspr.modelling.MultipleModelsStarter;
import qspr.modelling.SetStatistics;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.modelling.configurations.ExternalCondition;
import qspr.util.ClassCompressor;
import qspr.util.RWriter;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.useractions.EventFactory;
import com.eadmet.useractions.ModelCreationAction;
import com.eadmet.utils.NumericalValueStandardizer;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.mailer.Mailer;

public class MultipleModelsService 
{

	private static final int MAXIMUM_DESCRIPTORS_LENGTH = 60;

	private static final Logger logger = LogManager.getLogger(MultipleModelsService.class);

	private String getObjectMD5(Serializable object)
	{
		if(object == null) return "";
		return OCHEMUtils.getMD5(ClassCompressor.objectToByte(object));
	}

	private void delete(Model model, ModelData mData)
	{
		logger.info("Deleting model " + model.id);
		model.delete();
		mData.status = ModelData.STATUS.deleted;
	}

	private CalculationClient getClient()
	{
		CalculationClient cc = new CalculationClient("Model summary");
		cc.setDeepSleepTime(1);
		cc.setTolerateMetaserverDown();
		return cc;
	}

	public List<MultipleModelsData> getMultipleModelsData(Long setId, Boolean deleteInvalidModels) throws Exception
	{
		//desccache cachedesc deschash hashdesc

		MultipleModelsData mmDataCV = new MultipleModelsData().setValidation(QSPRConstants.CV, "Cross-Validation");
		MultipleModelsData mmDataBagging = new MultipleModelsData().setValidation(QSPRConstants.BAGGING, QSPRConstants.BAGGING);
		MultipleModelsData mmDataNoValidation = new MultipleModelsData().setValidation(QSPRConstants.NO_VALIDATION, "Other protocols");
		Basket trainingSet = mmDataCV.trainingSet = mmDataBagging.trainingSet = mmDataNoValidation.trainingSet = (Basket) Globals.session().get(Basket.class, setId);

		PendingTaskPeer.updateTaskStatuses(trainingSet.user, trainingSet, 200);
		List<Model> models = Repository.model.getModelsByTrainingSet(trainingSet);

		int processed = 0;
		long lastRestart = System.nanoTime(); // Last time a transaction was restarted.
		Set<Long> uniqueIDs = new HashSet<Long>();

		for (Model model : models)
			uniqueIDs.add(model.id);

		int modelsUnique = uniqueIDs.size();
		uniqueIDs.clear();

		List<Integer> setSizes = new ArrayList<Integer>();

		boolean excluded = false;

		for (Model model : models)
		{
			if (uniqueIDs.contains(model.id))
				continue;
			uniqueIDs.add(model.id);
			processed++;
			model = (Model) Globals.session().get(Model.class, model.id);
			model.recalculateModelSize();
			if (model.template.isDescriptorCalculationOnly())
				continue;
			Object conf;
			StandartizationOptions standartization = null;

			try {
				conf = model.attachment.getObject().configuration;
				standartization =  model.attachment.getObject().standartization;
			}catch(Exception e) {
				String message = "Obsolete template: cannot display model with template: " + model.template.name + " for model: " + model.publicId + " " + model.name;
				Mailer.notifyAdmins("Obsolete template detected", message);
				continue;
			}
			CDSConfiguration cdsConf = null;
			if (conf instanceof CDSConfiguration)
				cdsConf = (CDSConfiguration) conf;

			MultipleModelsData mmData;
			if (model.attachment.getObject().protocol.validationConfiguration == null)
				mmData = mmDataNoValidation;
			else
				mmData = (model.attachment.getObject().protocol.validationConfiguration instanceof BaggingConfiguration) ? mmDataBagging : mmDataCV;

			mmData.hasAnyModels = true;

			logger.info("Processing model " + processed + " out of " + modelsUnique + ": " + model.name +" id: "+model.id+ " taskId: "+model.taskId);

			String methodHash,methodName,descHash,descName;


			if (cdsConf != null)
			{
				descHash = getObjectMD5(cdsConf.descriptors);
				if (cdsConf.selection.getDescriptorsSize() > 0)
					descHash += getObjectMD5(cdsConf.selection); // Include the descriptors selection into the descriptor configuration hash

				if(!StructureOptimisationConfiguration.isDefaultConfiguration(cdsConf.optimisationConfiguration))
					descHash += cdsConf.optimisationConfiguration.getTaskType();

				if(cdsConf.conditions != null) {
					String s = "";
					for(ExternalCondition d: cdsConf.conditions)
						s += d.toString();
					descHash += getObjectMD5(s);
				}


				// Skipping large uploaded model
				CompressedObject<Object> uploadedModelData = cdsConf.modelConfiguration.uploadedModelData;
				if(uploadedModelData != null)
					cdsConf.modelConfiguration.uploadedModelData = null; 
				ModelAbstractConfiguration duplicatedConfig = cdsConf.modelConfiguration.getDeepCopy();
				cdsConf.modelConfiguration.uploadedModelData = uploadedModelData;
				// Skipping large uploaded model


				String  aug = "";

				if(duplicatedConfig instanceof NoDescriptors) { // required to avoid the use of augmentation a part of the method 
					NoDescriptors nodesc = (NoDescriptors)duplicatedConfig;
					aug = nodesc.augemenationString();
					nodesc.setAugmentations(1, 1, false);
				}

				String version = duplicatedConfig.versionOCHEM;
				if(version == null || version.startsWith("v.")){
					duplicatedConfig.setVersion("", true); // we do not need to distinguish different OCHEM versions unless this is manually set version
					version = "";
				}

				methodHash = getObjectMD5(duplicatedConfig);

				//desc_hash
				//hash_desc

				methodName = duplicatedConfig.getInformativeName()+ (version.length()>0?" (" + version+ ")":"");
				String stand = standartization == null || standartization.getDefault() == ChemInfEngine.CDK?"":" ("+ standartization.getDefault()+")";
				
				methodName += stand;
				methodHash += stand;
				descName = cdsConf.descriptors.types.toString().replaceAll("[\\[\\]]", "");

				if(descName.length() == 0) {
					descHash += model.attachment.getObject().standartization.desaltWith != null;
					descName = "SMILES";
					if(model.attachment.getObject().standartization.desaltWith == null)descName += " (mix)";
				}

				descName += aug;
				descHash += aug;
			}
			else
			{
				if (!QSPRConstants.CONSENSUS.equals(model.template.name) && !QSPRConstants.UPLOADED_MODEL.equals(model.template.name) )
					continue;
				ConsensusModelConfiguration cons = (ConsensusModelConfiguration) conf;
				descHash = cons.toString();
				descName = QSPRConstants.CONSENSUS + " ("+cons.individualModels.size()+" models)";
				methodHash = model.template.name + " " + cons.type;
				methodName = QSPRConstants.CONSENSUS+ ":" + cons.type;
			}


			if(descName.length() > MAXIMUM_DESCRIPTORS_LENGTH) 
				descName = descName.substring(0, MAXIMUM_DESCRIPTORS_LENGTH) + "...";

			if(cdsConf != null && cdsConf.selection != null && cdsConf.selection.useUFS)
				descName += " (ufs)";

			if(cdsConf != null && cdsConf.optimisationConfiguration != null) {
				String name = cdsConf.optimisationConfiguration.getTaskType().toLowerCase();
				if(!name.equals("molconvertor"))
					descName += " 3D:" + name;
			}

			if(conf instanceof ProvidedConditions && ((ProvidedConditions) conf).hasConditions())
				descName += " + Cond.";

			if(cdsConf != null && cdsConf.hasMixtures())
				if(cdsConf.descriptors.mixtures == MixturesProcessing.FRACTION)
					descName += " (mix)";
				else
					descName += " + Mix:" + cdsConf.descriptors.mixtures;

			descHash += descName;

			if (model.modelMappings.get(0).statisticsRecalculated == null)
				continue;

			if (!mmData.descriptorHashes.contains(descHash))
			{
				mmData.descriptorHashes.add(descHash);
				mmData.descriptorCodes.add(descName);
			}

			String trainingSetHash = model.getTrainingSetHash();

			methodHash += trainingSetHash;

			if(model.isStratifyValidationConfiguration()) methodHash += "str";

			if(model.attachment.getObject().protocol.validationConfiguration instanceof ValidationConfiguration) {
				ValidationConfiguration valid = (ValidationConfiguration)model.attachment.getObject().protocol.validationConfiguration;
				methodHash += valid.ensembleSize;
				methodHash += valid.getInformativeName();
				methodName += valid.getInformativeName();
			}

			if (!mmData.trainingSetHashes.contains(trainingSetHash))
				mmData.trainingSetHashes.add(trainingSetHash);

			int tsNum = mmData.trainingSetHashes.indexOf(trainingSetHash);
			if (!mmData.methodHashes.contains(methodHash))
			{
				mmData.methodHashes.add(methodHash);

				if (tsNum > 0)
					methodName += " (tr. set. "+(tsNum + 1)+")";
				String methodCo = methodName;
				int k = 1;
				while (mmData.methodCodes.contains(methodCo))
					methodCo = methodName + "(" + (++k) + ")";
				mmData.methodCodes.add(methodCo);
			}

			ModelData mData = new ModelData();
			mData.setModel(model);

			List<PendingTask> pTasks = PendingTask.getByModel(model, TaskType.MODEL_TRAINING);
			PendingTask pTask = null;
			if (!pTasks.isEmpty())
			{
				pTask = pTasks.get(0);
				mData.pendingTaskId = pTask.id;
			}

			boolean special = model.template.name.equals(QSPRConstants.CONSENSUS) || model.template.name.equals(QSPRConstants.UPLOADED_MODEL);

			if (!special && !model.isStatisticsCalculated) // model is not yet finished / fetched
				if (pTasks.isEmpty())
					mData.status = ModelData.STATUS.invalid;
				else
				{
					mData.setStatus(pTask.status);
					mData.detailedStatus = pTask.getDetailedStatus();
					mData.pendingTaskId = pTask.id;

					if ("ready".equals(pTask.status))try
					{
						ModelProcessor processor = ModelFactory.getProcessor(model.template);
						processor.model = model;
						processor.onTaskReceived(getClient().getTask(pTask.taskId));
						processor.saveModel();
					}
					catch (Exception e)
					{
						logger.info("Failed to save model " + model.name);
						e.printStackTrace();
					}
				}

			if (model.isStatisticsCalculated || special) // model is finished
				for(int j = 0; j < model.modelMappings.size(); j++)
				{
					ModelMapping mapping = model.modelMappings.get(j); //

					ModelStatistics ms = (ModelStatistics)mapping.statisticsRecalculated.getObject();
					ms.setBootstrapReplicas(QSPRConstants.NO_REPLICAS);
					ms.recalculateStatistics(mapping);
					mData.status = ModelData.STATUS.ready;
					if (model.taskId == null)
						mData.status = model.published ? ModelData.STATUS.published : ModelData.STATUS.saved;
					if (model.modelMappings.isEmpty())
					{
						mData.status = ModelData.STATUS.invalid;
						mData.detailedStatus = "No model mapping found! The model is invalid or has been damaged";
					}
					else
					{
						List<SetStatistics> sets = order(ms.sets, setSizes);

						boolean first = j == 0;
						for (SetStatistics ss : sets){
							mData.addStatistics(ss, model, first);
							first = false;
						}

						if(sets.get(sets.size()-1).setId.equals(QSPRConstants.EXCLUDED)){
							SetStatistics s = new SetStatistics(sets.get(0),sets.get(sets.size()-1));
							s.bootstrapReplicas = QSPRConstants.NO_REPLICAS;
							s.recalculateStatistics(mapping);
							s.setId = QSPRConstants.WHOLE;
							mData.addStatistics(s, model, false);
							excluded = true;
						}else{
							mData.addStatistics(sets.get(0), model, false); // to avoid recalculation of statistics
							mData.stats.get(mData.stats.size()-1).id = QSPRConstants.WHOLE; 
						}
					}
				}

			if(mmData.modelData == null || mmData.modelData.length <= mmData.methodHashes.indexOf(methodHash))
				mmData.addMethod();

			if(mmData.modelData[0].length <= mmData.descriptorHashes.indexOf(descHash))
				mmData.addDescriptors();

			if (mmData.modelData[mmData.methodHashes.indexOf(methodHash)][mmData.descriptorHashes.indexOf(descHash)] != null)
				mmData.modelData[mmData.methodHashes.indexOf(methodHash)][mmData.descriptorHashes.indexOf(descHash)].addAnotherModel(mData);
			else
				mmData.modelData[mmData.methodHashes.indexOf(methodHash)][mmData.descriptorHashes.indexOf(descHash)] = mData;

			if (ModelData.STATUS.invalid == mData.status && deleteInvalidModels)
				delete(model, mData);

			if (System.nanoTime() - lastRestart > 5 * 1E9)
			{
				logger.info("Restarting transactions after "+ NumericalValueStandardizer.getSignificantDigits((System.nanoTime() - lastRestart) / 1E9) + " seconds of activity");
				lastRestart = System.nanoTime();
				Globals.restartAllTransactions(true);
				logger.info("done");
			}
		}

		mmDataCV.sortDescriptors();
		mmDataBagging.sortDescriptors();
		mmDataNoValidation.sortDescriptors();

		if(!excluded){
			mmDataCV.cleanWhole();
			mmDataBagging.cleanWhole();
			mmDataNoValidation.cleanWhole();
		}

		mmDataCV.anonymyseValidationSets();
		mmDataBagging.anonymyseValidationSets();
		mmDataNoValidation.anonymyseValidationSets();

		mmDataCV.addAverage();
		mmDataBagging.addAverage();
		mmDataNoValidation.addAverage();

		List<MultipleModelsData> data = new ArrayList<MultipleModelsData>();
		data.add(mmDataCV);
		data.add(mmDataBagging);
		data.add(mmDataNoValidation);
		return data;
	}

	void printSizes(MultipleModelsData mmData ) {

		if(mmData.modelData == null) 
			System.out.println("Length: 0");
		else {
			System.out.println("Length: "+ mmData.modelData.length);
			for(ModelData[] m : mmData.modelData) 
				System.out.println(" -- "+ m.length);
		}		
	}

	/**
	 * Order validation sets to be always in the same order
	 * @param sets
	 * @param setSizes
	 * @return
	 */

	private List<SetStatistics> order(List<SetStatistics> sets,
			List<Integer> setSizes) {

		List<SetStatistics> ordered = new ArrayList<SetStatistics>();

		SetStatistics last = null;

		for(SetStatistics s:sets)
			if(QSPRConstants.TRAINING.equals(s.setId))
				ordered.add(s);
			else
				if(QSPRConstants.EXCLUDED.equals(s.setId))
					last = s;
				else
					if(!setSizes.contains(s.points.size()))setSizes.add(s.points.size());

		for(Integer i:setSizes)
			for(SetStatistics s:sets)
				if(s.setId.startsWith(QSPRConstants.VALIDATION) && s.points.size() == i.intValue())
					ordered.add(s);

		if(last != null)
			ordered.add(last);

		return ordered;
	}

	public ExportableModel createCrossOverModel(Long methodId, Long descrId)
	{
		Model meth = Repository.model.getById(methodId);
		Model descr = Repository.model.getById(descrId);

		ExportableModel eModel = ExportableModel.create(descr); // copy of descriptors

		CDSConfiguration cds = (CDSConfiguration) eModel.attachment.configuration;
		ModelAbstractConfiguration modelcfg = ((CDSConfiguration) meth.attachment.getObject().configuration).modelConfiguration;

		if(!cds.isCompatibleDescriptorsAndMethod(modelcfg))return null;

		// copy from Method
		eModel.method = meth.template.name;
		cds.modelConfiguration = modelcfg; // copy all required
		eModel.attachment.protocol = meth.attachment.getObject().protocol;
		eModel.attachment.datahandling = meth.attachment.getObject().datahandling;
		eModel.attachment.standartization = meth.attachment.getObject().standartization;

		eModel.attachment.standartization.synchronise();

		//However deSalt - which defines mixture comes from descriptors
		if(descr.attachment.getObject().standartization.desaltWith == null)
			eModel.attachment.standartization.desaltWith = null;
		else
			eModel.attachment.standartization.desaltWith = eModel.attachment.standartization.getDefault(); // keeping the same but with different method

		CDSConfiguration descConf = (CDSConfiguration) descr.attachment.getObject().configuration;  // original desc configuration which was overwritten
		cds.selection = descConf.selection.getDeepCopy(); // selection comes from descriptors

		if(cds.modelConfiguration instanceof NoDescriptors){ // also augmentation
			cds.modelConfiguration = cds.modelConfiguration.getDeepCopy(); // otherwise the same modelConfiguration template could be used again  
			NoDescriptors conf = (NoDescriptors) descConf.modelConfiguration;
			((NoDescriptors)cds.modelConfiguration).setAugmentations(conf.getAugementTraining(), conf.getAugmentApply(), conf.getBalanceData());
		}

		eModel.attachment.standartization.synchronise();
		return eModel;
	}

	/*

	 	public ExportableModel createCrossOverModel(Long methodModelId, Long descriptorModelId)
	{
		Model modelMethod = Repository.model.getById(methodModelId);
		Model modelDescriptors = Repository.model.getById(descriptorModelId);

		ExportableModel eModel = ExportableModel.create(modelDescriptors);
		eModel.method = modelMethod.template.name;
		CDSConfiguration eModelOutConf = (CDSConfiguration) eModel.attachment.configuration;
		CDSConfiguration eModelInDesc = (CDSConfiguration) modelDescriptors.attachment.getObject().configuration;

		if(eModelOutConf.modelConfiguration instanceof NoDescriptors != eModelInDesc.modelConfiguration instanceof NoDescriptors) return null; // both should be of the same type

		eModelOutConf.modelConfiguration = ((CDSConfiguration) modelMethod.attachment.getObject().configuration).modelConfiguration;
		eModel.attachment.protocol = modelMethod.attachment.getObject().protocol;
		if(eModel.attachment.protocol.validationConfiguration != null && eModel.attachment.protocol.validationConfiguration.byMaximumComponent())
			eModel.attachment.standartization.desaltWith = modelMethod.attachment.getObject().standartization.desaltWith;
		eModel.attachment.datahandling = modelMethod.attachment.getObject().datahandling;
		eModelOutConf.selection = eModelInDesc.selection.getCopy();
		eModelOutConf.selection.descriptors = eModelInDesc.selection.descriptors;

		if(eModelOutConf.modelConfiguration instanceof NoDescriptors){
			eModelOutConf.modelConfiguration = eModelOutConf.modelConfiguration.getDeepCopy(); // otherwise the same modelConfiguration template is always used 
			NoDescriptors model = (NoDescriptors)eModelOutConf.modelConfiguration;
			model.setAugmentations(1, 1, false); // first remove any augmentation
			if(eModelInDesc.modelConfiguration instanceof NoDescriptors){ // it should be like this!!!
				NoDescriptors conf = (NoDescriptors) eModelInDesc.modelConfiguration;
				model.setAugmentations(conf.getAugementTraining(), conf.getAugmentApply(), conf.getBalanceData());
			}
		}

		return eModel;
	}

	 */

	public MultipleModelsStarter createMultipleCrossOverModels(List<Long> methodModelIds, List<Long> descriptorModelIds, boolean savedOnly) throws Exception
	{
		if (methodModelIds.size() != descriptorModelIds.size())
			throw new UserFriendlyException("The count of method templates and descriptor templates does not match");

		int nModels = methodModelIds.size();

		MultipleModelsStarter starter = new MultipleModelsStarter();
		starter.trainingSetId = Repository.model.getById(Long.valueOf(methodModelIds.get(0))).trainingSet.id;
		for (int i = 0; i < nModels; i++)
		{
			if(savedOnly) {
				Model m = Repository.model.getById(descriptorModelIds.get(i));
				List<PendingTask> pTasks = PendingTask.getByModel(m, TaskType.MODEL_TRAINING);
				if(pTasks.size() != 0)continue;
			}

			ExportableModel mod = createCrossOverModel(methodModelIds.get(i), descriptorModelIds.get(i));
			if(mod != null)starter.eModels.add(mod);

			Model m = Repository.model.getById(methodModelIds.get(i));
			for (ModelMapping mm : m.modelMappings)
				starter.units.put(mm.property, mm.unit);
		}

		starter.skipDublicates = false;

		starter.start();
		return starter;
	}

	@SuppressWarnings("unchecked")
	public List<ModelConfigurationTemplate> getModelTemplates()
	{
		return Globals.session().createCriteria(ModelConfigurationTemplate.class).add(Restrictions.ne("type", TemplateType.MODEL))
				.add(Restrictions.or(Restrictions.eq("introducer", Globals.userSession().user), Restrictions.eq("isPublic", true)))
				.list();
	}

	public MultipleModelsStarter createMultipleModels(CrossModelGenerator generator)
	{
		MultipleModelsStarter modelStarter = new MultipleModelsStarter();

		modelStarter.trainingSetId = generator.trainingSetId;
		modelStarter.validationSetId = generator.validationSetIds;

		for (int i=0; i<generator.propertyIds.size(); i++)
			modelStarter.units.put(Property.getById(generator.propertyIds.get(i)), Unit.getById(generator.unitIds.get(i)));

		modelStarter.skipDublicates = generator.skipDublicates;
		modelStarter.eModels = generator.generateModelTemplates();
		modelStarter.start();

		EventFactory.document("Multiple models creation", new ModelCreationAction(Basket.getBasket(Globals.userSession(), modelStarter.trainingSetId), modelStarter.eModels.size()), null);
		return modelStarter;
	}

	public void exportAsR(MultipleModelsData[] models, BufferedOutputStream os) throws IOException
	{
		RWriter writer = new RWriter(os);

		for (MultipleModelsData mmData : models) 
		{
			writer.variablePrefix = "models$" + mmData.validation + "$";
			for (int i = 0; i < mmData.methodCodes.size(); i++) {
				for (int k = 0; k < mmData.descriptorCodes.size(); k++)
				{
					if (mmData.modelData[i][k] != null && mmData.modelData[i][k].stats != null && mmData.modelData[i][k].stats.size() > 0)
					{
						switch(mmData.getModelType()) {
						case classification:
							writer.addValue("accuracy", mmData.modelData[i][k].stats.get(0).accuracy); 
							writer.addValue("balancedAccuracy", mmData.modelData[i][k].stats.get(0).balancedAccuracy);
							writer.addValue("auc", mmData.modelData[i][k].stats.get(0).auc); 
							break;
						case mixed:
							writer.addValue("rmse", mmData.modelData[i][k].stats.get(0).rmse); // TODO: Consider if we have more than 1 model in a cell // TODO: Consider validation sets
							writer.addValue("r2", mmData.modelData[i][k].stats.get(0).r2);
							writer.addValue("q2", mmData.modelData[i][k].stats.get(0).q2);
							writer.addValue("mae", mmData.modelData[i][k].stats.get(0).mae);
							writer.addValue("accuracy", mmData.modelData[i][k].stats.get(0).accuracy); 
							writer.addValue("balancedAccuracy", mmData.modelData[i][k].stats.get(0).balancedAccuracy);
							writer.addValue("auc", mmData.modelData[i][k].stats.get(0).auc); 
							break;
						case regression:
							writer.addValue("rmse", mmData.modelData[i][k].stats.get(0).rmse); // TODO: Consider if we have more than 1 model in a cell // TODO: Consider validation sets
							writer.addValue("r2", mmData.modelData[i][k].stats.get(0).r2);
							writer.addValue("q2", mmData.modelData[i][k].stats.get(0).q2);
							writer.addValue("mae", mmData.modelData[i][k].stats.get(0).mae);
							break;
						}
					}
					else
					{
						writer.addNA("rmse");
						writer.addNA("r2");
						writer.addNA("q2");
						writer.addNA("mae");

						writer.addNA("accuracy");
						writer.addNA("balancedAccuracy");
						writer.addNA("auc");
					}
				}
			}

			for (int i = 0; i < mmData.methodCodes.size(); i++)
				writer.addString("methodNames",  mmData.methodCodes.get(i));
			for (int i = 0; i < mmData.descriptorCodes.size(); i++)
				writer.addString("descriptorNames",  mmData.descriptorCodes.get(i));

			writer.write();

			writer.writeDirect("reorder = function(m) m[order(rownames(m)),]");
			String tmp = "models$val$__ = matrix(models$val$__, ncol = length(models$val$methodNames), nrow = length(models$val$descriptorNames), dimnames = list(models$val$descriptorNames, models$val$methodNames))\n";
			tmp += "models$val$__ = reorder(models$val$__)";
			tmp = tmp.replaceAll("val", mmData.validation);
			writer.writeDirect(tmp.replaceAll("__", "rmse"));
			writer.writeDirect(tmp.replaceAll("__", "r2"));
			writer.writeDirect(tmp.replaceAll("__", "q2"));
			writer.writeDirect(tmp.replaceAll("__", "mae"));
		}

		writer.close();
	}

	public void exportAsXls(MultipleModelsData[] mmDatas, BufferedOutputStream os) throws IOException
	{
		String[] metrics = mmDatas[0].getModelType() == ModelType.regression ? new String[]{"rmse", "r2", "q2", "mae", "average", "modelIDs"} :
			mmDatas[0].getModelType() == ModelType.classification ?
					new String[]{"auc","accuracy", "balancedAccuracy", "modelIDs"} :
						new String[]{"rmse", "r2", "q2", "mae", "auc", "accuracy", "balancedAccuracy", "average", "modelIDs"}
		;

		// Enumerate all basket IDs in this CM summary
		LinkedHashSet<String> basketIDs = new LinkedHashSet<String>();
		for (MultipleModelsData mmData : mmDatas)
			for (int i = 0; i < mmData.descriptorCodes.size(); i++)
				for (int k = 0; k < mmData.methodCodes.size(); k++)
					if (mmData.modelData[k][i] != null)
						for (int j = 0; j < mmData.modelData[k][i].stats.size(); j++)
							basketIDs.add(mmData.modelData[k][i].stats.get(j).getBasketSetNameId());

		// Create a sheet per validation type + basket + metrics
		HSSFWorkbook workbook = new HSSFWorkbook();
		for (MultipleModelsData mmData : mmDatas) 
		{
			int bNum = 0;
			for (String basketID : basketIDs)
			{
				Long id = Long.valueOf(basketID.split("_")[0]);
				String suffix = basketID.split("_")[1];
				switch(suffix) {
				case "training": suffix = "tr"; break;
				case "validation": suffix = "val" + bNum++; break;
				case "excluded": suffix = "excl"; break;
				}

				Basket basket = (Basket) Globals.session().get(Basket.class, id);

				ModelData m = null;	String sep = null;

				for (int i = 0; i < mmData.descriptorCodes.size(); i++) 
					for (int k = 0; k < mmData.methodCodes.size(); k++)if(mmData.modelData[i][k]!=null && mmData.modelData[i][k].stats !=null )
					{
						m=mmData.modelData[i][k];
						i=mmData.descriptorCodes.size();
						break;
					}

				if(m != null)
					for(SetStats ss: m.stats) 
						if(ss != null && ss.getBasketSetNameId().equals(basketID))
							sep = sep == null?"":sep+";";

				if(sep==null)sep="";

				for (String metric : metrics) 
				{
					HSSFSheet sheet = workbook.createSheet(mmData.validation + "-" + suffix + "-" + metric);
					boolean sheetEmpty = true;
					sheet.createRow(0).createCell(0).setCellValue(basket != null ? basket.name : "Unknown basket");
					HSSFRow headerRow = sheet.createRow(1);
					for (int k = 0; k < mmData.methodCodes.size(); k++)
						headerRow.createCell(k + 1).setCellValue(mmData.methodCodes.get(k)+sep);
					for (int i = 0; i < mmData.descriptorCodes.size(); i++) 
					{
						HSSFRow row = sheet.createRow(i + 2);
						row.createCell(0).setCellValue(mmData.descriptorCodes.get(i).replaceAll(";", "+"));
						for (int k = 0; k < mmData.methodCodes.size(); k++)
						{
							if (mmData.modelData[k][i] != null)
							{
								if ("modelIDs".equals(metric))
									row.createCell(k + 1).setCellValue(mmData.modelData[k][i].publicModelId);
								else
								{
									ModelData mm = mmData.modelData[k][i];
									String all = null;
									for(SetStats ss: mm.stats) 
										if(ss.getBasketSetNameId().equals(basketID)) {
											String val = null;
											switch(metric) {
											case "accuracy": val = trimParentheses(ss.accuracy); break;
											case "balancedAccuracy": val = trimParentheses(ss.balancedAccuracy); break;
											case "rmse": val = trimParentheses(ss.rmse); break;
											case "r2": val = trimParentheses(ss.r2); break;
											case "q2": val = trimParentheses(ss.q2); break;
											case "mae": val = trimParentheses(ss.mae); break;
											case "average": val = trimParentheses(ss.average); break;
											case "auc": val = trimParentheses(ss.auc); break;
											default: throw new UserFriendlyException("Non defined metric: "+ metric);
											}
											if(val != null)
												all = all != null ? all + ";" + val : val;
										}

									if (all != null)
									{
										sheetEmpty = false;
										if(all.contains(";"))
											row.createCell(k + 1).setCellValue(all);
										else
											row.createCell(k + 1).setCellValue(Double.parseDouble(all));
									}
								}
							}
						}
					}

					if (sheetEmpty)
						workbook.removeSheetAt(workbook.getSheetIndex(sheet));

				}
			}
		}
		workbook.write(os);
	}

	String trimParentheses(String v) {
		if( v == null)return null;
		return v.replace("(", "").replace(")", "");
	}

}
