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

package qspr.controllers;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.ThreadScope;
import qspr.dao.Repository;
import qspr.entities.Alert;
import qspr.entities.Model;
import qspr.entities.ModelConfigurationTemplate;
import qspr.entities.ModelConfigurationTemplate.TemplateType;
import qspr.export.ExportableModel;
import qspr.frontend.MultipleModelsData;
import qspr.frontend.WebModel;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.StandartizationOptions;
import qspr.metaserver.configurations.DescriptorsConfiguration.MixturesProcessing;
import qspr.metaserver.configurations.ModelAbstractConfiguration.ScalingType;
import qspr.metaserver.configurations.ValidationConfiguration.MixtureValidation;
import qspr.modelling.CrossModelGenerator;
import qspr.modelling.DataPreprocessingParser;
import qspr.modelling.MultipleModelsStarter;
import qspr.modelling.configurators.ConditionConfigurator;
import qspr.modelling.configurators.DescriptorsConfigurator;
import qspr.modelling.configurators.ModelConfigurationTemplateFactory;
import qspr.util.AccessChecker;
import qspr.util.ClassCompressor;

import com.eadmet.business.MultipleModelsService;
import com.eadmet.exceptions.UserFriendlyException;

/**
 * Controller for the Comprehensive Modeling utility
 * Comprehensive Modeling allows to analyze multiple models for the same endpoint and training set
 * @author midnighter
 *
 */
@Controller
public class MultipleModelsController extends ControllerWrapper 
{
	MultipleModelsService service = new MultipleModelsService();

	public MultipleModelsController()
	{
		sessionRequired = true;
	}
	/**
	 * Show the tabular model overview screen
	 * @throws Exception
	 */
	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) throws Exception
	{

		ThreadScope.get().considerPredicatesInStatistics = assertParam("consider-predicates");

		List<MultipleModelsData> data = service.getMultipleModelsData(getLongParam("set"), assertParam("deleteInvalidModels"));

		Globals.setSessionAttribute(SessionVariable.MULTIPLE_MODELS_REPORT, data.toArray(new MultipleModelsData[]{}));

		WebModel wm = new WebModel();
		for (MultipleModelsData multipleModelsData : data) 
			if (multipleModelsData.hasAnyModels)
				wm.addObject(multipleModelsData);

		if(ThreadScope.get().considerPredicatesInStatistics) // indicate that we did not considered predicates
			wm.addParam("consider-predicates", "consider-predicates");

		return wm
				.setTemplate("model/multiple-models")
				.getModelAndView();
	}

	/**
	 * Create a model from an "empty" cell in the table (click "+")
	 */
	public ModelAndView createCrossOverModel(HttpServletRequest request, HttpServletResponse response)
	{
		ExportableModel eModel = service.createCrossOverModel(getLongParam("method-from"), getLongParam("descriptors-from"));
		if(eModel == null)throw new UserFriendlyException("Configurations are not compatible");
		Globals.setSessionAttribute(SessionVariable.IMPORTED_TEMPLATE, eModel);
		return redirect(getRedirectCheckConsensus(eModel));
	}

	/**
	 * Create multiple "missing" models in empty table cells
	 * @throws Exception 
	 */
	public ModelAndView createMultipleCrossOverModels(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ClassCompressor.classLoader = MultipleModelsController.class.getClassLoader();//?

		String[] strMethodFrom = getParam("method-from").split(",");
		String[] strDescriptorFrom = getParam("descriptors-from").split(",");
		String saved = getParam("saved");

		List<Long> methodFrom = new ArrayList<Long>();
		List<Long> descriptorFrom = new ArrayList<Long>();

		for (String strMethod : strMethodFrom)
			methodFrom.add(Long.valueOf(strMethod));

		for (String strDescriptor : strDescriptorFrom)
			descriptorFrom.add(Long.valueOf(strDescriptor));

		service.createMultipleCrossOverModels(methodFrom, descriptorFrom,"true".equals(saved));

		return new WebModel().getModelAndView();
	}

	/**
	 * Show the comprehensive modeling screen
	 * @return
	 */
	public ModelAndView create(HttpServletRequest request, HttpServletResponse response)
	{
		StandartizationOptions options = new StandartizationOptions();
		options.setDefaults();
		return new WebModel(options).addObjects(service.getModelTemplates()).setTemplate("model/create-multiple-models").getModelAndView();
	}

	/**
	 * Show the comprehensive modeling screen
	 * @return
	 */
	public ModelAndView add(HttpServletRequest request, HttpServletResponse response)
	{
		return new WebModel().addObjects(service.getModelTemplates()).setTemplate("model/add-model-template").getModelAndView();
	}


	/**
	 * Submit the dialog and initiate creation of multiple models based on the selected templates
	 * @throws IOException 
	 */
	public ModelAndView createSubmit(HttpServletRequest request, HttpServletResponse response) throws IOException
	{
		CrossModelGenerator generator = new CrossModelGenerator();
		String[] sids = getParam("ids").split(",");

		generator.skipDublicates = assertParam("skipDublicates");

		for (String sid : sids)
			generator.templateIds.add(Integer.valueOf(sid));

		generator.trainingSetId = getLongParam("trainingsetid");
		generator.validationSetIds = ModelController.getValidationSetIDs(request);

		if(generator.validationSetIds != null)
			for(Long v: generator.validationSetIds)
				if( (long)v == (long)generator.trainingSetId)
					throw new UserFriendlyException("The training set should not be used as a validation set");

		int i = 0;
		do
		{
			if (!assertParam("property" + i + "-id"))
				break;
			generator.propertyIds.add(getLongParam("property" + i + "-id"));
			generator.unitIds.add(getLongParam("unit" + i));
			i++;
		} while (true);

		DataPreprocessingParser.parseStandartizationUI(request, generator.standartization, generator.dataHandling, null);

		if (assertParam("override-normalisation"))
		{
			ScalingType[] scaling = DataPreprocessingParser.parseNormalisationUI(request);
			generator.scalingTypeX = scaling[0];
			generator.scalingTypeY = scaling[1];
		}

		generator.optimisationConfiguration = DescriptorsConfigurator.configureStructureOptimisation(request);

		if(generator.standartization.desaltWith == null) { // mixture is activated
			generator.mixturesProcessing = MixturesProcessing.fromString(getParam("mixtures"));
			if(getParam("mixtureValidation") != null) {
				generator.mixtureValidation = MixtureValidation.valueOf(getParam("mixtureValidation").toUpperCase());
				if(generator.mixtureValidation == MixtureValidation.MIXTURE) generator.mixtureValidation = null;
			}
		}
		generator.saveModels = assertParam("save-models");
		generator.saveModels = generator.saveModels ? null : false;

		generator.additionaDescriptors = getParamTrim("additionalDescriptorsMerge");
		if(generator.additionaDescriptors != null && generator.additionaDescriptors.length()> 0 && !generator.additionaDescriptors.contains(":"))
			generator.additionaDescriptors=Globals.getUsername()+":"+generator.additionaDescriptors;

		generator.singleLearning = assertParam("singleLearningEnforce");
		generator.featureNet = assertParam("featureNet");
		generator.sanitize = assertParam("sanitize");
		generator.crs = assertParam("crs");
		if(generator.additionaDescriptors != null && generator.additionaDescriptors.length() == 0)generator.additionaDescriptors = null;
		generator.externalDescriptors = ConditionConfigurator.addConditions(request);
		generator.implicitValues = getParamTrim("implicitValues");
		generator.taxonomy = getParamTrim("taxonomyValues");
		if(generator.implicitValues != null && generator.implicitValues.length() ==0)generator.implicitValues = null;
		generator.versionOCHEM = getVersionInfo();

		MultipleModelsStarter modelStarter = service.createMultipleModels(generator);

		Globals.setSessionAttribute(SessionVariable.MULTIPLE_MODELS_STARTER, modelStarter);
		return new WebModel().getModelAndView();
	}



	public static String getDescriptorsCode(DescriptorsConfiguration descConf)
	{
		return descConf.types.toString();
	}

	public ModelAndView getModelXml(HttpServletRequest request, HttpServletResponse response) throws IOException, ClassNotFoundException, Exception
	{
		Model model = Repository.model.getById(getLongParam("model"));
		return new WebModel(new Alert(model.configurationXml)).getModelAndView();
	}


	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws IOException, ClassNotFoundException, Exception
	{
		String[] ids = getParam("modelId").split(",");
		for (String id : ids) 
		{
			Model model = (Model) Globals.session().get(Model.class, Long.valueOf(id));
			if (model == null)
				continue; // The model has already been deleted 

			if(model.published)
				throw new UserFriendlyException("The model is published and can't be modified.");

			if (!model.getPrivileges(request).canEdit) 
				throw new UserFriendlyException("You do not have sufficient privileges for this model" + AccessChecker.explainUser(Repository.user.getBySessionId(model.session.id)));

			switch(getParam("action")) {
			case "delete": model.delete(); break;
			case "terminate": terminateTasks(model); break;
			case "save": 
				model.storePendingModel();
				model.updateDescription();
				break;
			}
		}

		return new WebModel().getModelAndView();
	}

	public ModelAndView exportReportAsR(HttpServletRequest request, HttpServletResponse response) throws IOException
	{
		MultipleModelsData[] mmDatas = (MultipleModelsData[]) Globals.getSessionAttribute(SessionVariable.MULTIPLE_MODELS_REPORT);
		if (mmDatas == null)
			throw new UserFriendlyException("Report data is not available. Probably you have been logged out. Please, reload the page with the model summary");

		response.setContentType("text/plain");
		response.setHeader("Content-Disposition", "attachment; filename=models_report_"+mmDatas[0].trainingSet.name.replaceAll("\\s+", "_")+".R");
		BufferedOutputStream os = new BufferedOutputStream(response.getOutputStream());
		service.exportAsR(mmDatas, os);
		os.flush();
		return null;
	}

	public ModelAndView exportReportAsXls(HttpServletRequest request, HttpServletResponse response) throws IOException
	{
		MultipleModelsData[] mmDatas = (MultipleModelsData[]) Globals.getSessionAttribute(SessionVariable.MULTIPLE_MODELS_REPORT);
		if (mmDatas == null)
			throw new UserFriendlyException("Report data is not available. Probably you have been logged out. Please, reload the page with the model summary");

		response.setContentType("text/plain");
		response.setHeader("Content-Disposition", "attachment; filename=models_report_"+mmDatas[0].trainingSet.name.replaceAll("\\s+", "_")+".xls");
		BufferedOutputStream os = new BufferedOutputStream(response.getOutputStream());
		service.exportAsXls(mmDatas, os);
		os.flush();
		return null;
	}

	public ModelAndView addCustomTemplate(HttpServletRequest request, HttpServletResponse response) throws IOException
	{
		ModelConfigurationTemplate mcTemplate = ModelConfigurationTemplateFactory.create(Repository.model.getById(getLongParam("model")), TemplateType.valueOf(getParam("type")));
		mcTemplate.isPublic = false;
		mcTemplate.session = Globals.userSession();
		mcTemplate.introducer = Globals.userSession().user;
		mcTemplate.name = "Created from " + Repository.model.getById(getLongParam("model")).name;
		mcTemplate.updateHash();

		mcTemplate.checkConflicts();
		Globals.session().save(mcTemplate);

		return new WebModel().setTemplate("model/create-multiple-models").getModelAndView();
	}

	@SuppressWarnings("unchecked")
	private void terminateTasks(Model model) throws IOException, ClassNotFoundException
	{
		CalculationClient client = getClient();
		List<Integer> taskIDs = Globals.session().createQuery("select taskId from PendingTask where model=:model and (status='init' or status='assigned')").setParameter("model", model).list();
		for (Integer taskID : taskIDs) 
			client.killTask(taskID);
	}

	private CalculationClient getClient()
	{
		CalculationClient cc = new CalculationClient("Model summary");
		cc.setDeepSleepTime(1);
		cc.setTolerateMetaserverDown();
		return cc;
	}

}
