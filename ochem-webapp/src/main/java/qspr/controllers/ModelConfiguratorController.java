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

import java.io.File;
import java.lang.reflect.InvocationTargetException;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Hibernate;
import org.hibernate.criterion.Order;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.dao.Repository;
import qspr.entities.Alert;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.ModelConfigurationTemplate;
import qspr.entities.ModelMapping;
import qspr.entities.ModelTemplate;
import qspr.entities.Property;
import qspr.entities.Unit;
import qspr.export.ExportableModel;
import qspr.frontend.WebModel;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.configurations.BaggingConfiguration;
import qspr.metaserver.configurations.CrossValidationConfiguration;
import qspr.metaserver.configurations.ValidationConfiguration;
import qspr.metaserver.configurations.ValidationConfiguration.MixtureValidation;
import qspr.modelling.configurators.BasicModelConfigurator;
import qspr.util.BasicRecordMapper;
import qspr.workflow.utils.QSPRConstants;

@Controller
public class ModelConfiguratorController extends ControllerWrapper
{
	private static transient final Logger logger = LogManager.getLogger(ModelConfiguratorController.class);

	public ModelConfiguratorController()
	{
		sessionRequired = true;
	}

	public ModelAndView configure(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		BasicModelConfigurator configurator;
		Model model = null;
		if (assertParam("model"))
		{
			model = (Model) Globals.session().get(Model.class, getLongParam("model"));
			configurator = BasicModelConfigurator.getConfigurator(model.template);
			configurator.model = model;
			Globals.setSessionAttribute(SessionVariable.MODEL_CONFIGURATOR, configurator);
			Globals.setSessionAttribute(SessionVariable.MODEL, model);
		}
		else
			configurator = getConfigurator();
		if (assertParam("page"))
			configurator.currentPage = getParam("page");
		else
		{
			if (configurator.upload)
				return redirect("modelconfigurator/configure.do?page=" + configurator.currentPage);
			else
				return redirect("modelconfigurator/configure.do?page=" + configurator.currentPage + "&upload=1");
		}

		configurator.versionOCHEM = getVersionInfo();

		WebModel web = (WebModel) configurator.getClass().getMethod(configurator.currentPage).invoke(configurator);

		model = (Model)Globals.getSessionAttribute(SessionVariable.MODEL);

		web.addParam("use-mixtures-validation", model != null? model.canBeMixtures():"true");
		web.addParam("use-mixtures-descriptors", model.attachment.getObject().standartization.desaltWith == null?"true":"false");

		return web.addParam("page", configurator.currentPage)
				.addParam("upload", "1").getModelAndView();
	}

	public ModelAndView importTemplate(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		return new WebModel().setTemplate("modeller/import").getModelAndView();
	}

	public ModelAndView importTemplateSubmit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		File f = Globals.getUploadedFile();

		ExportableModel eModel = (ExportableModel) Globals.jaxbContext.createUnmarshaller().unmarshal(f);
		Globals.setSessionAttribute(SessionVariable.IMPORTED_TEMPLATE, eModel);

		return WebModel.redirect(getRedirectCheckConsensus(eModel)).getModelAndView();
	}

	public ModelAndView createFromTemplate(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ExportableModel template = null;
		if (assertParam("model"))
			template = ExportableModel.create(Repository.model.getById(getLongParam("model")));
		if (assertParam("template"))
			template = (ExportableModel) ModelConfigurationTemplate.getById(getIntParam("template")).configuration.getObject();
		if (assertParam("trainingSet"))
		{
			Basket basket = Basket.getBasket(Globals.userSession(), getLongParam("trainingSet"));
			template.trainingSetId = basket.id;
		}
		Globals.setSessionAttribute(SessionVariable.IMPORTED_TEMPLATE, template);
		return WebModel.redirect(getRedirectCheckConsensus(template)).getModelAndView();
	}

	public ModelAndView status(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		BasicModelConfigurator configurator = getConfigurator();
		return new WebModel(new Alert(configurator.getStatus())).getModelAndView();
	}

	public ModelAndView kill(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		BasicModelConfigurator configurator = getConfigurator();
		new CalculationClient("Web-Interface").killTask(configurator.model.taskId);
		return new WebModel().getModelAndView();
	}

	public ModelAndView submit(HttpServletRequest request, HttpServletResponse response) throws Throwable
	{
		BasicModelConfigurator configurator = getConfigurator();

		if (assertParam("page"))
			configurator.currentPage = getParam("page");
		try
		{
			configurator.getClass().getMethod(configurator.currentPage + "Submit", HttpServletRequest.class).invoke(configurator, request);
		} catch (InvocationTargetException e)
		{
			throw e.getCause();
		}

		if (configurator.upload)
			return redirect("modelconfigurator/configure.do?upload=1&page=" + configurator.currentPage);
		else
			return redirect("modelconfigurator/configure.do?page=" + configurator.currentPage);
	}

	public ModelAndView choose(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ExportableModel template = (ExportableModel) Globals.getSessionAttribute(SessionVariable.IMPORTED_TEMPLATE);
		return new WebModel(template).setList(Globals.session().createCriteria(ModelTemplate.class).addOrder(Order.asc("name")).list())
				.setTemplate(assertParam("descriptors") ? "modeller/select-training-set-descriptors" : "modeller/select-training-set").getModelAndView();
	}

	public ModelAndView choosesubmit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Model model = new Model();

		if (assertParam("skip-configuration"))
			model = ((ExportableModel) Globals.getSessionAttribute(SessionVariable.IMPORTED_TEMPLATE)).createModel();

		BasicRecordMapper brm;
		model.session = Globals.userSession();
		if (assertParam("trainingsetid"))
			model.trainingSet = Basket.getBasket(Globals.userSession(), getLongParam("trainingsetid"));

		if(model.trainingSet == null || model.trainingSet.id == null)
			throw new UserFriendlyException("No training set was provided or available: select a training set for your model.");
		
		long entries = Repository.basket.countEntries(model.trainingSet.id);

		Repository.user.checkEligibility(entries, QSPRConstants.MODEL_BONUS);

		List<Long> validationSetIDs = ModelController.getValidationSetIDs(request);
		model.cleanValidationSets();
		if (validationSetIDs != null && !validationSetIDs.isEmpty())
			for (Long valSetID : validationSetIDs)
			{
				model.addValidationSet(Basket.getBasket(Globals.userSession(), valSetID));
			}

		// Configure Multi Learning Classes ("model mappings")
		brm = new BasicRecordMapper(model.trainingSet);
		List<Property> propList = model.trainingSet.getProperty();
		int count = 0;
		model.modelMappings.clear();
		for (Property property : propList)
		{
			Hibernate.initialize(property.options);
			ModelMapping modelMapping = new ModelMapping();
			modelMapping.property = property;

			if (assertParam("unit" + count))
				modelMapping.unit = (Unit) Globals.session().get(Unit.class, getLongParam("unit" + count));
			else
				modelMapping.unit = property.defaultUnit;
			modelMapping._class = brm.getClass(modelMapping.property);
			model.addModelMapping(modelMapping);
			count++;
		}

		BasicModelConfigurator configurator;
		String firstStep;
		if (assertParam("skip-configuration"))
		{
			// We will use the template and skip the detailed model configuration
			logger.info("Using the template configuration:\n" + model.attachment.getObject());
			configurator = BasicModelConfigurator.getConfigurator(model.template);
			configurator.model = model;
			firstStep = "start";
		}
		else
		{

			Long templateID = getLongParam("template");
			logger.info("getLongParam(\"template\") returns: " + templateID);
			if (templateID == null) {
				for (Object param : request().getParameterMap().keySet()) {
					String param2 = getParam((String)param);
					if (param2.length() >0) {
						logger.info("parameter: " + param + "\t\thas value: " + param2);
					}
				}
			}
			model.template = (ModelTemplate) Globals.session().get(ModelTemplate.class, templateID);
			configurator = BasicModelConfigurator.getConfigurator(model.template);
			configurator.importedTemplate = (ExportableModel) Globals.getSessionAttribute(SessionVariable.IMPORTED_TEMPLATE);

			configurator.setModel(model);

			if (assertParam("upload"))
				configurator.upload = true;

			model.attachment.getObject().protocol.validationConfiguration = null; // be default no validation

			ValidationConfiguration config = null;

			if (QSPRConstants.BAGGING.equalsIgnoreCase(getParam(QSPRConstants.VALIDATION))){
					config = new BaggingConfiguration();
			}
			if (QSPRConstants.CV.equalsIgnoreCase(getParam(QSPRConstants.VALIDATION)))
				config = new CrossValidationConfiguration();

			model.attachment.getObject().protocol.validationConfiguration = config;

			if (config != null)
			{ // we have validation protocol
				try {
				
				System.out.println("valid: " + getParam(QSPRConstants.VALIDATION));
				System.out.println("valid-num: " + getParam(getParam(QSPRConstants.VALIDATION) + "-ensemble"));

				if (getParam("record-stratified") != null)
					config.mixtureValidation = MixtureValidation.RECORD;
				//if (getParam("record-desault") != null) {
				//	if(config.recordValidated != null && config.recordValidated)
				//		throw new UserFriendlyException("Validations by records and by maximum component are not comaptible");
				//	config.deSaulted = true;
				//}

				config.ensembleSize = Integer.parseInt(getParam(getParam(QSPRConstants.VALIDATION) + "-ensemble"));
				//if (assertParam(getParam(QSPRConstants.VALIDATION) + "-mixture-validated"))
				//	config.mixtureValidation = true; //TO CHANGE
				if (getParam(getParam(QSPRConstants.VALIDATION) + "-stratified") != null)
					config.validationType = BaggingConfiguration.STRATIFIED;
				if (getParam("bagging-instances") != null && config instanceof BaggingConfiguration)
					((BaggingConfiguration) config).numberInstances = Integer.parseInt(getParam("bagging-instances"));
				if (getParam("store-individual-predictions") != null && config instanceof BaggingConfiguration)
					((BaggingConfiguration) config).keepIndividualPredictions = true;
				}catch(Exception e) {
					throw new UserFriendlyException("A temporal failure (presumably a problem with the cache of the browser). Please, try again.");
				}
			}

			model.isCompatibleModelAndDescriptors();
			
			firstStep = configurator.dataPreprocessingStep ? "dataPreprocessing" : configurator.firstConfigurationStep;
		}

		Globals.setSessionAttribute(SessionVariable.MODEL_CONFIGURATOR, configurator);
		Globals.setSessionAttribute(SessionVariable.MODEL, model);

		//model.checkInconsistencies();

		if (configurator.upload)
			return redirect("modelconfigurator/configure.do?page=" + firstStep + "&upload=1");
		else
			return redirect("modelconfigurator/configure.do?page=" + firstStep);
	}

	protected ModelAndView handleRequestInternal(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		try
		{
			return super.handleRequestInternal(request, response);
		} catch (InstantiationException e)
		{
			return redirect("modelconfigurator/choose.do");
		}
	}

	private BasicModelConfigurator getConfigurator() throws Exception
	{
		BasicModelConfigurator object = (BasicModelConfigurator) Globals.getSessionAttribute(SessionVariable.MODEL_CONFIGURATOR);
		if (object == null)
			throw new InstantiationException();
		return object;
	}
}