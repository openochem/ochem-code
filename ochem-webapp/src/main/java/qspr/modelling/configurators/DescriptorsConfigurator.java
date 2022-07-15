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

package qspr.modelling.configurators;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.servlet.http.HttpServletRequest;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.HibernateException;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.business.toxalert.AlertsFilter;
import qspr.business.toxalert.ScreeningProcessor;
import qspr.dao.Repository;
import qspr.entities.Alert;
import qspr.entities.CalculatedDescriptor;
import qspr.entities.CalculatedDescriptor.CalculatedDescriptorStatus;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.Property;
import qspr.entities.SubstructureAlert;
import qspr.entities.User;
import qspr.frontend.AvailableAlertsData;
import qspr.frontend.WebModel;
import qspr.interfaces.ProvidedConditions;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.configurations.OpenBabelConfiguration;
import qspr.metaserver.configurations.StandardNoDescriptorsConfiguration;
import qspr.metaserver.configurations.BalloonConfiguration;
import qspr.metaserver.configurations.CorinaConfiguration;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorType;
import qspr.metaserver.configurations.DescriptorsConfiguration.MixturesProcessing;
import qspr.metaserver.configurations.DescriptorsExpValuesConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.metaserver.configurations.ModelAbstractConfiguration.ScalingType;
import qspr.metaserver.configurations.NoDescriptors;
import qspr.metaserver.configurations.ObgenConfiguration;
import qspr.metaserver.configurations.DescriptorsApplyModelConfiguration;
import qspr.metaserver.configurations.DescriptorsStructuralAlertsConfiguration;
import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.workflow.utils.QSPRConstants;
import qspr.modelling.CDSModelProcessor;
import qspr.modelling.DataPreprocessingParser;
import qspr.modelling.ModelFactory;
import qspr.modelling.ModelProcessor;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.modelling.configurations.ExternalCondition;
import qspr.util.BasicRecordMapper;
import qspr.util.WrapperThread;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.NodeConfiguration;
import qspr.workflow.datatypes.StringList;
import qspr.workflow.datatypes.WorkflowConfiguration;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.descriptorcache.DescriptorConfigEntry;
import com.eadmet.descriptorcache.DescriptorsRepository;
import com.eadmet.descriptorcache.DescriptorsRepositoryFactory;
import com.eadmet.exceptions.UserFriendlyException;

public class DescriptorsConfigurator extends BasicModelConfigurator
{
	private static transient final Logger logger = LogManager.getLogger(DescriptorsConfigurator.class);

	public static final int DEFAULT_3D_OPTIMIZATION = 1; // CORINA with 3D
	public int recommendedOptimisation = 0;
	public String stepAfterDescriptors = startStep;
	List<String> descriptorOverrides = new ArrayList<String>();
	public List<CalculatedDescriptor> manualDescriptorList = new ArrayList<CalculatedDescriptor>();

	public DescriptorsConfigurator()
	{
		this.firstConfigurationStep = "descriptorBlocks";
	}

	public WebModel calculateDescriptors()
	{
		asynkTask = new WrapperThread()
		{

			public void setStatus(String status)
			{
				if (!status.equals(this.getStatus()))
				{
					super.setStatus(status);
					logger.info(status);
				}
			}

			public void wrapped() throws Exception
			{
				CDSConfiguration configuration = (CDSConfiguration) model.attachment.getObject().configuration;

				model.mergeTrainingAndValidationSet();

				ModelProcessor mp = ModelFactory.getProcessor(model.template);
				mp.model = model;
				mp.mapper = new BasicRecordMapper(model.trainingSet);

				WorkflowNodeData wndInput = (WorkflowNodeData) mp.getDataForTeacher();
				WorkflowConfiguration wc = new WorkflowConfiguration("CDS")
						.addNodeConfiguration(new NodeConfiguration("descriptors-processor", configuration.descriptors))
						.addNodeConfiguration(new NodeConfiguration("descriptors-selector", configuration.selection));

				CDSModelProcessor.addMoleculeProcessingConfiguration(wc.nodesConfiguration, model.attachment.getObject().standartization, configuration.optimisationConfiguration);

				Task task = new Task(QSPRConstants.Workflow, wc, wndInput,false);

				task.setPriority(TaskPriority.HIGH);
				CalculationClient client = new CalculationClient("Web-interface/"+Globals.getUsername(),Globals.getUsername())
				{
					public void setStatus(String taskStatus)
					{
						super.setStatus(taskStatus);
						if (taskStatus != null)
							if (!taskStatus.equals(asynkTask.getStatus()))
								asynkTask.setStatus(taskStatus);
					}
				};

				task = client.calculateTask(task);
				task.check();
				setStatus("Parsing client response");
				DataTable dtDescriptors = ((WorkflowNodeData) task.getResult()).ports.get(0);
				processCalculatedDescriptors(dtDescriptors); // required to set descriptors that should be preselected

				setStatus("Finished");
			}
		};
		asynkTask.userSession = Globals.userSession();

		asynkTask.start();

		return new WebModel().setTemplate("modeller/configurators/calculate-descriptors");
	}

	public void calculateDescriptorsSubmit(HttpServletRequest request)
	{
		asynkTask = null;
		this.currentPage = "selectDescriptors";
	}

	public WebModel descriptorBlocks()
	{
		boolean useAtomic = false;

		if (model != null)
			for (ModelMapping mapping : model.modelMappings)
				if (mapping.property.getShortName().contains("pka"))
				{
					useAtomic = true;
					break;
				}

		WebModel wm = new WebModel(model).addParam("UseAtomic", useAtomic ? "true" : "false").setTemplate("modeller/configurators/descriptor-configuration-step");

		if (importedTemplate != null)
			wm.addObject(importedTemplate);

		// Mixture Related Code
		wm.addParam("use-mixtures-descriptors", model == null || model.attachment.getObject().standartization.desaltWith != null?"false":"true");
		// Mixture Related Code End

		// Get the list of supported tasks
		try
		{
			Set<String> taskTypes = new CalculationClient("DescriptorsModelConfigurator").getSupportedTaskTypes();

			for (String taskType : taskTypes)
				wm.addParam("available-task", taskType);

		} catch (Exception e)
		{
			logger.error("Could not load supported task types", e);
			wm.addParam("availability-unknown", "true");
		}

		List<Model> featuredModels = Repository.model.getFeaturedModels();
		wm.addObjects(featuredModels);

		try
		{
			User user = Globals.userSession().user;
			DescriptorsRepository dRepository = DescriptorsRepositoryFactory.getRepository();
			if(user != null)
				wm.addObjects(dRepository.getConfigurations(user.login));

			wm.addObjects(dRepository.getConfigurations(QSPRConstants.PUBLISHER)); // public data
		}
		catch (Exception e)
		{
			logger.error("Could not request descriptor configuration from the storage", e);
		}

		// Info required for the alerts filter
		wm.addObject(new AvailableAlertsData());

		return wm;
	}

	public void descriptorBlocksSubmit(HttpServletRequest request) throws Exception
	{
		CDSConfiguration configuration = (CDSConfiguration) model.attachment.getObject().configuration;
		configureDescriptors(request, configuration);

		// Fail-safe check that at least one descriptor was selected
		if (!configuration.hasDescriptors() && !configuration.hasConditions())
		{
			if(configuration.modelConfiguration instanceof NoDescriptors) 
				this.currentPage = this.stepAfterDescriptors;
			else
				this.currentPage = "descriptorBlocks";
			return;
		}

		this.currentPage = "structureOptimisation";

		// / skip optimization if it's not required... ugly dublication of code
		// - but oh well... //NoS 28.06.10
		if (recommendedOptimisation == 0)
		{
			//			configuration.optimisationConfiguration = StructureOptimisationConfiguration.getDefaultConfiguration();
			if (upload)
				this.currentPage = "uploadFile";
			else
				this.currentPage = "descriptorFilters";
		}
	}

	void addFeatureNets(HttpServletRequest request, CDSConfiguration configuration)
	{
		PredictionScenario scenario = null;
		if (request.getParameter("predictionScenario") != null)
			scenario = PredictionScenario.valueOf(request.getParameter("predictionScenario"));
		String[] models_id = request.getParameterValues("model-id");
		for (int i = 0; i < models_id.length; i++)
		{
			Model innerModel = (Model) Globals.session().get(Model.class, Long.valueOf(models_id[i]));

			String columnTitle = "";
			if (innerModel.modelMappings.size() > 1)
				columnTitle = "Multi-Classes model";
			else
				columnTitle = innerModel.modelMappings.get(0).property.getName() + "(" + innerModel.name + ")";

			DescriptorsApplyModelConfiguration predictionConfig = new DescriptorsApplyModelConfiguration();
			predictionConfig.modelId = innerModel.publicId;
			predictionConfig.scenario = scenario;
			configuration.descriptors.addDescriptorType(DescriptorsConfiguration.ApplyModel + "_" + innerModel.publicId, predictionConfig, columnTitle, null,true);

			ProvidedConditions conf = (ProvidedConditions)innerModel.attachment.getObject().configuration;

			if(conf.hasConditions()) { 
				configuration.conditions = ExternalCondition.checkAndMerge(
						configuration.conditions,
						conf.getConditions(), innerModel.publicId);

				model.attachment.getObject().optionsMapping = innerModel.attachment.getObject().optionsMapping;
			}
		}
	}

	DescriptorsStructuralAlertsConfiguration addStructuralAlerts(HttpServletRequest request) throws HibernateException, Exception
	{
		AlertsFilter filter = new AlertsFilter();
		filter.parseUI(request);
		@SuppressWarnings("unchecked")
		List<SubstructureAlert> alerts = filter.filterCriteria().list();
		logger.info("Using " + alerts.size() + " alerts");
		DescriptorsStructuralAlertsConfiguration conf = new DescriptorsStructuralAlertsConfiguration();
		conf.compactMode = false;
		conf.alertPatterns = ScreeningProcessor.getSMARTs(alerts);
		return conf;
	}


	public CDSConfiguration configureDescriptors(HttpServletRequest request, CDSConfiguration configuration) throws Exception
	{
		DescriptorsConfiguration descriptors = configuration.descriptors; 

		descriptors.allowMerge = request.getParameter("allow_merge") != null? true : null;
		descriptors.forceUpdateDescriptorCache = request.getParameter("force_cache") != null ? true : null;
		descriptors.mixtures = MixturesProcessing.fromString(request.getParameter("mixtures"));

		// add "standard" descriptors
		descriptors.setConfiguration(request);

		// first stored ones are added
		addStoredDescriptors(request, descriptors);

		// Few methods requires an access to Properties/Conditions are implemented here
		if (request.getParameterValues("model-id") != null)
			addFeatureNets(request, configuration);

		if (request.getParameter(DescriptorsConfiguration.StructuralAlerts) != null)
			descriptors.addDescriptorType(DescriptorsConfiguration.StructuralAlerts, addStructuralAlerts(request));

		if (request.getParameter(DescriptorsConfiguration.ExpValues) != null)
			descriptors.addDescriptorType(DescriptorsConfiguration.ExpValues, addExpValues(request));

		if(configuration.modelConfiguration == null || configuration.modelConfiguration.isSupportConditions())
			configuration.conditions = ConditionConfigurator.addConditions(request);

		recommendedOptimisation = descriptors.requires3D()? DEFAULT_3D_OPTIMIZATION : 0;

		return configuration;
	}

	private void addStoredDescriptors(HttpServletRequest request, DescriptorsConfiguration descriptors) throws Exception {
		// Stored descriptors
		for (Object key : request.getParameterMap().keySet())
		{
			DescriptorsRepository dRepository = DescriptorsRepositoryFactory.getRepository();
			String param = key.toString();
			if (param.startsWith("Stored-"))
			{
				String objectID = param.substring(7);
				DescriptorConfigEntry dConfig = dRepository.getDescriptorConfigById(objectID);

				DescriptorsAbstractConfiguration descConfiguration = null;
				if (dConfig.description != null && !dConfig.description.equals(""))
					descConfiguration = (DescriptorsAbstractConfiguration) Globals.jaxbContext.createUnmarshaller().unmarshal(new StringReader(dConfig.description));

				DescriptorType dType = descriptors.addDescriptorType(dConfig.type, descConfiguration,null,null,true);
				dType.markUncachedAsErrors = true;
			}
		}		
	}

	private DescriptorsAbstractConfiguration addExpValues(HttpServletRequest request) {

		DescriptorsExpValuesConfiguration eConf = new DescriptorsExpValuesConfiguration();

		eConf.multipleExpValues = request.getParameter("exp-values-multiple");
		eConf.basketId = Long.valueOf(request.getParameter("exp-values-basket-id"));
		String[] propIds = request.getParameterValues("exp-values-property-id");
		for (String propId : propIds)
		{
			Property property = Property.getById(Long.valueOf(propId));
			eConf.properties.add(property.getName());
		}

		return eConf;
	}


	public WebModel descriptorFilters()
	{
		WebModel wm = new WebModel(getDefaultTemplate()).setTemplate("modeller/configurators/descriptor-filters").addParam("multipart", "true");
		if (model != null && model.template != null)
			if (Arrays.asList(new String[]{QSPRConstants.ASNN, QSPRConstants.MLRA, QSPRConstants.KNN}).contains(model.template.name))
				wm.addParam("no-normalisation", "true");
		descriptorOverrides.clear();
		return wm;
	}

	private List<String> parseDescriptorListFile(File f) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(f));
		try
		{
			String line;
			List<String> descriptorOverrides = new ArrayList<String>();
			while ((line = reader.readLine()) != null)
				descriptorOverrides.add(line.trim());
			return descriptorOverrides;
		} finally
		{
			reader.close();
		}
	}

	public void descriptorFiltersSubmit(HttpServletRequest request)
	{
		CDSConfiguration configuration = (CDSConfiguration) model.attachment.getObject().configuration;
		HttpServletRequest mp = ThreadScope.get().localRequest;

		try
		{
			File fileWithDescriptorList = Globals.getUploadedFile("file-with-descriptor-list");
			if (fileWithDescriptorList != null)
				descriptorOverrides = parseDescriptorListFile(fileWithDescriptorList);
		} catch (Exception e)
		{
			e.printStackTrace();
			request.getSession().setAttribute("desc-list-error", e.getMessage());
			this.currentPage = "descriptorFilters";
			return;
		}

		if (mp.getParameter("unique-values") != null)
			configuration.selection.numDifferentValues = Integer.parseInt(mp.getParameter("unique-values-threshold"));
		else
			configuration.selection.numDifferentValues = 0;

		if (mp.getParameter("decorrelation") != null)
			configuration.selection.correlationThreshold = Double.parseDouble(mp.getParameter("decorrelation-threshold"));
		else
			configuration.selection.correlationThreshold = 0;

		if (mp.getParameter("maximum") != null)
		{
			long v = (int) Long.parseLong(mp.getParameter("maximum-threshold"));
			if (v > Integer.MAX_VALUE)
				v = Integer.MAX_VALUE;
			configuration.selection.maximumValueThreshold = (int) v;
		}
		else
			configuration.selection.maximumValueThreshold = Integer.MAX_VALUE;

		if (mp.getParameter("std") != null)
			configuration.selection.stdThreshold = Double.parseDouble(mp.getParameter("std-threshold"));
		else
			configuration.selection.stdThreshold = 0;

		if (mp.getParameter("ufs") != null)
			configuration.selection.useUFS = true;

		if (((CDSConfiguration) model.attachment.getObject().configuration).modelConfiguration != null)
		{
			ModelAbstractConfiguration conf = ((CDSConfiguration) model.attachment.getObject().configuration).modelConfiguration;

			ScalingType[] scaling = DataPreprocessingParser.parseNormalisationUI(request);
			if(conf.isSupportDescriptors())conf.scaleTypeX = scaling[0];
			conf.scaleTypeY = scaling[1];
		}

		if (mp.getParameter("manual") != null)
			this.currentPage = "calculateDescriptors";
		else
			this.currentPage = stepAfterDescriptors;
	}

	public WebModel previewFile()
	{
		return new WebModel().setTemplate("modeller/configurators/preview-file").setObject(model.preview);
	}

	public void previewFileSubmit(HttpServletRequest request)
	{
		this.currentPage = "calculateDescriptors";
	}

	private void processCalculatedDescriptors(DataTable dtDescriptors) throws IOException
	{
		long id = 0;
		manualDescriptorList.clear();
		Map<String, CalculatedDescriptor> local = new HashMap<String, CalculatedDescriptor>();

		for (String column : dtDescriptors.getColumns()) 
		{
			CalculatedDescriptor descriptor = new CalculatedDescriptor();
			descriptor.id = --id;
			if (descriptorOverrides.size() == 0) //Select all descriptors only if there are no overrides available
				descriptor.selected = true;

			descriptor.name = column;			
			manualDescriptorList.add(descriptor);
			local.put(column, descriptor);
		}

		for (String override : descriptorOverrides) 
		{
			if (local.containsKey(override))
			{
				if (local.get(override).status.equals(CalculatedDescriptorStatus.NORMAL))
				{
					local.get(override).status = CalculatedDescriptorStatus.EXPECTED_FOUND;
					local.get(override).selected = true;
				} //Else  - no change
			}
			else
			{
				CalculatedDescriptor descriptor = new CalculatedDescriptor();
				descriptor.id = --id;
				descriptor.selected = true;
				descriptor.name = override;
				descriptor.status = CalculatedDescriptorStatus.EXPECTED_NOT_FOUND;
				manualDescriptorList.add(descriptor);
				local.put(override, descriptor);
			}
		}


	}

	/**
	 * Required to provide model upload from a file
	 * 
	 * @param f
	 * @throws Exception
	 */
	public void processModelUpload(File f) throws Exception
	{
	}

	public WebModel selectDescriptors()
	{
		WebModel wm = new WebModel(importedTemplate).setTemplate("modeller/configurators/select-descriptors");
		if (descriptorOverrides.size() > 0)
			wm.addParam("show-selected", "true");
		return wm;
	}

	public void selectDescriptorsSubmit(HttpServletRequest request) throws IOException
	{
		StringList sl = new StringList();

		for (CalculatedDescriptor d : manualDescriptorList)
			if (d.selected)
				sl.values.add(d.name);

		CDSConfiguration configuration = (CDSConfiguration) model.attachment.getObject().configuration;
		configuration.selection.storeDescriptorAsStrings(sl.values);

		if (configuration.selection.getDescriptorsSize() > 0)
			this.currentPage = stepAfterDescriptors;
	}

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		if (model.attachment.getObject().configuration == null)
			model.attachment.getObject().configuration = new CDSConfiguration();
	}

	public WebModel structureOptimisation()
	{
		// Hardcode recommendation to optimize here
		WebModel wm = new WebModel(model).addParam("recommendation", String.valueOf(recommendedOptimisation)).setTemplate(
				"modeller/configurators/structure-optimisation");
		if (importedTemplate != null)
			wm.addObject(importedTemplate);
		return wm;
	}

	static public StructureOptimisationConfiguration configureStructureOptimisation(HttpServletRequest request)
	{

		if(request.getParameter("optimisation") == null) return new BalloonConfiguration(); 
		
		Integer optimisationOption = Integer.valueOf(request.getParameter("optimisation"));

		switch (optimisationOption)
		{
		case 1:
			return new CorinaConfiguration();
		case 4:
			return new OpenBabelConfiguration();
		case 5:
			return new ObgenConfiguration();
		case 6:
			return new BalloonConfiguration();
		case 0:
			return null;
		default:
			throw new UserFriendlyException("structure optimisation is not determined for: " + optimisationOption);
		}
	}

	public void structureOptimisationSubmit(HttpServletRequest request)
	{
		CDSConfiguration configuration = (CDSConfiguration) model.attachment.getObject().configuration;

		configuration.optimisationConfiguration = configureStructureOptimisation(request);

		if (upload)
			this.currentPage = "uploadFile";
		else
			this.currentPage = "descriptorFilters";
	}

	private String message;

	public WebModel uploadFile()
	{
		return new WebModel(new Alert(message)).setTemplate("modeller/configurators/upload-file").addParam("multipart", "true");
	}

	/**
	 * Uploading a model
	 * 
	 * @param request
	 */

	public void uploadFileSubmit(HttpServletRequest request) // DescriptorsModelConfigurator.java
	{
		CDSConfiguration configuration = (CDSConfiguration) model.attachment.getObject().configuration;
		File f = null;

		try
		{
			f = Globals.getUploadedFile();
		} catch (Exception e)
		{
			request.getSession().setAttribute("desc-list-error", e.getMessage());
			this.currentPage = "uploadFile";
			return;
		}

		try
		{
			processModelUpload(f);
			message = "";
		} catch (Exception e)
		{
			message = e.getMessage();
			request.getSession().setAttribute("desc-list-error", e.getMessage());
			this.currentPage = "uploadFile";
			return;
		}

		configuration.selection.setNoFiltering(); // in order to skip some filtering of descriptors

		// configuration.selection.descriptors = new StringList();
		this.currentPage = "previewFile";
	}

	static void addAugmentation(NoDescriptors conf, HttpServletRequest request) {
		Integer validation = null, training = request.getParameter("augmentation") == null ? null :Integer.valueOf(request.getParameter("augmentation"));
		validation = request.getParameter("augmentApplySet") == null ? null :Integer.valueOf(request.getParameter("augmentApplySet"));
		if(training != null && training > QSPRConstants.MAX_AUGMENTATION) training = QSPRConstants.MAX_AUGMENTATION;
		if(validation != null && validation > QSPRConstants.MAX_AUGMENTATION) validation = QSPRConstants.MAX_AUGMENTATION;
		conf.setAugmentations(training, validation,request.getParameter("balance") != null);
	}

	static void addStandard(StandardNoDescriptorsConfiguration conf, HttpServletRequest request, int maxepochs) {
		conf.shuffle =  request.getParameter("shuffle") != null;
		conf.nepochs = Integer.valueOf(request.getParameter("nepochs"));
		if(request.getParameter("early") != null) conf.early = Double.valueOf(request.getParameter("early"));
		conf.chirality = request.getParameter("chirality") != null; // option is to ignore it
		if(request.getParameter("sanitize") != null)
			conf.setSanitize();

		conf.nepochs = conf.nepochs < 5? 5: conf.nepochs > maxepochs ? maxepochs:conf.nepochs;
		conf.early = conf.early == null || conf.early < 0.05? 0.05: conf.early >1?1: conf.early;

		conf.activationFunction = request.getParameter("activation_function") == null ? null : 
			StandardNoDescriptorsConfiguration.ACTIVATE.valueOf(request.getParameter("activation_function"));

		if(conf.activationFunction == StandardNoDescriptorsConfiguration.ACTIVATE.RELU) conf.activationFunction = null;

		addAugmentation(conf, request);
	}

}
