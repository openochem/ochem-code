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

import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;

import org.hibernate.Criteria;
import org.hibernate.Hibernate;
import org.hibernate.criterion.ProjectionList;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.SessionVariable;
import qspr.business.ModelPeer;
import qspr.dao.Repository;
import qspr.entities.Attachment;
import qspr.entities.Attachment.AttachmentType;
import qspr.entities.AttachmentSource;
import qspr.entities.Basket;
import qspr.entities.DataHandlingOptions;
import qspr.entities.Model;
import qspr.entities.ModelAttachment;
import qspr.entities.ModelMapping;
import qspr.entities.ModelMicroAttachment;
import qspr.entities.ModelProtocol;
import qspr.entities.ModelTemplate;
import qspr.entities.PendingTask;
import qspr.entities.Property;
import qspr.entities.PendingTask.TaskType;
import qspr.export.ExportableModel;
import qspr.frontend.WebModel;
import qspr.metaserver.configurations.LabelWeighting;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.configurations.StandartizationOptions;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.modelling.DataPreprocessingParser;
import qspr.modelling.ModelFactory;
import qspr.modelling.ModelProcessor;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.util.WrapperThread;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.useractions.EventFactory;
import com.eadmet.useractions.ModelCreationAction;
import com.eadmet.useractions.ModelCreationAction.ModelAction;

public class BasicModelConfigurator
{
	public Model model;
	public ExportableModel importedTemplate;
	public String modelName;
	public String currentPage;
	public WrapperThread asynkTask;
	public boolean upload = false;
	public String startStep = "start";
	public String firstConfigurationStep = startStep;
	public boolean dataPreprocessingStep = true;
	public String preferredServer = OCHEMConfiguration.defaultPreferredServer;
	public String versionOCHEM;

	@SuppressWarnings("rawtypes")
	public static BasicModelConfigurator getConfigurator(ModelTemplate template)
	{
		if (template.name.equals("ASNN"))
			return new ANNConfigurator();
		try
		{
			Class configuratorClass = Class.forName("qspr.modelling.configurators." + template.name.replace("NEW", "").replace("-", "") + "Configurator");
			return (BasicModelConfigurator) configuratorClass.newInstance();
		} catch (Exception e)
		{
			return new BasicModelConfigurator();
		}
	}

	ExportableModel getDefaultTemplate(){
		if(importedTemplate != null) return importedTemplate;
		if(model != null) {
			ExportableModel eModel = ExportableModel.create(model);
			if(eModel.attachment.datahandling == null || model.id == null) eModel.attachment.datahandling = new DataHandlingOptions();
			if(eModel.attachment.standartization == null || model.id == null) eModel.attachment.standartization = new StandartizationOptions(true); // First model to create
			return eModel;
		}
		return null;
	}

	public BasicModelConfigurator()
	{
		this.currentPage = startStep;
	}

	public void setModel(Model model)
	{
		this.model = model;
		if (model.attachment == null)
			model.attachment = new Attachment<ModelAttachment>(new ModelAttachment(), AttachmentType.MARSHALABLE, AttachmentSource.Model);
		if (model.microattachment == null)
			model.microattachment = new Attachment<ModelMicroAttachment>(new ModelMicroAttachment(), AttachmentType.MARSHALABLE, AttachmentSource.Model);
		for (ModelMapping modelMapping : model.modelMappings)
			Hibernate.initialize(modelMapping.property.options);

		if (isWeightingRequired(model)) // So far only ANN will support label weighting
			startStep = "labelWeighting";
	}

	private boolean isWeightingRequired(Model model)
	{
		if (model.template == null) return false;

		if (!getDefaultWeighting(model.trainingSet).hasMultipleOutputs()) return false;

		return model.template.is_Support_Multilearning;

	}

	public WebModel start()
	{
		if (!model.template.isDescriptorCalculationOnly())
			Globals.setSessionAttribute(SessionVariable.IMPORTED_TEMPLATE, ExportableModel.create(model)); // remember the settings and use them next time
		model.name = ModelPeer.getModelName(model, true);
		return new WebModel(model).setTemplate("modeller/configurators/final-step");
	}

	public void startSubmit(HttpServletRequest request)
	{
		if (request.getParameter("discard") != null)
		{
			currentPage = "discard";
			return;
		}

		model.name = request.getParameter("name");
		model.defaultTaskPriority = Integer.valueOf(request.getParameter("priority"));

		// A "dedicated" task with a preferred server
		if (request.getParameter("preferred-server") != null && !request.getParameter("preferred-server").equals(""))
		{
			// If we have a dedicated task, assign a maximum priority
			model.defaultTaskPriority = TaskPriority.EXTRA_HIGH;
			preferredServer = request.getParameter("preferred-server").trim();
		}

		if (model.attachment.getObject().configuration instanceof CDSConfiguration)
			if (((CDSConfiguration) model.attachment.getObject().configuration).modelConfiguration != null)
			{
				ModelAbstractConfiguration conf = ((CDSConfiguration) model.attachment.getObject().configuration).modelConfiguration;
				conf.setVersion(versionOCHEM, false);
				if (request.getParameter("save-models") != null)
					conf.saveModels = null;
				else
					conf.saveModels = false;
			}

		currentPage = "teach";
		EventFactory.document("Model creation", new ModelCreationAction(model, ModelAction.START), null);
	}

	public WebModel teach()
	{
		if (asynkTask == null)
		{
			ModelProcessor processor = ModelFactory.getProcessor(model.template);
			processor.model = model;
			processor.preferredServer = preferredServer;
			processor.start();
			asynkTask = processor;
		}

		return new WebModel().setTemplate("modeller/configurators/teacher");
	}

	public void teachSubmit(HttpServletRequest request)
	{
		asynkTask = null;
		currentPage = "save";
	}

	public WebModel save() throws Exception
	{
		if ((model.taskId == null || model.taskId == 0) && model.id != null)
			model = (Model) Globals.session().get(Model.class, model.id); // Get a fresh copy from DB

		model.fetchCalculatedModel();

		if (model.name == null)
			model.name = "";
		model = (Model) Globals.session().merge(model);

		return new WebModel(model).setTemplate("modeller/savemodel");
	}

	public WebModel discard() throws Exception
	{
		return new WebModel().setTemplate("modeller/discard");

	}

	public void saveSubmit(HttpServletRequest request) throws Exception
	{
		model = (Model) Globals.session().get(Model.class, model.id);
		if (request.getParameter("discard") != null)
		{
			List<PendingTask> pTasks = PendingTask.getByModel(model, TaskType.MODEL_TRAINING);
			for (PendingTask pendingTask : pTasks)
				Globals.session().delete(pendingTask);
			if (!model.recalculation){
				model.delete();
			}
			currentPage = "discard";
			return;
		}

		if (!"".equals(request.getParameter("modelname")))
		{
			model.name = request.getParameter("modelname");
			//model.published = request.getParameter("published") != null;
			model.storePendingModel();
		}
		else
		{
			throw new UserFriendlyException("Empty model name provided");
		}

		EventFactory.document("Model save", new ModelCreationAction(model, ModelAction.SAVE));

		currentPage = "finished";
	}

	public WebModel later()
	{
		model = (Model) Globals.session().get(Model.class, model.id);

		if (asynkTask instanceof ModelProcessor)
			((ModelProcessor) asynkTask).exitAfterPosting();
		else
			asynkTask.interrupt();

		asynkTask = null;
		System.gc();
		return new WebModel(model).setTemplate("modeller/fetch-later");
	}

	public WebModel finished()
	{
		model = (Model) Globals.session().get(Model.class, model.id);
		return new WebModel(model).setTemplate("modeller/finished");
	}

	public String getStatus()
	{
		return asynkTask.getStatus();
	}

	@SuppressWarnings("unchecked")
	public WebModel dataPreprocessing()
	{
		Long numIntervals = 0L;
		Long numAppequals = 0L;
		Long numGreaterless = 0L;

		ProjectionList projList = Projections.projectionList();
		projList.add(Projections.groupProperty("pred.shortName"));
		projList.add(Projections.alias(Projections.countDistinct("exp.id"), "cnt"));

		List<Long> basketIdList = new ArrayList<Long>();
		basketIdList.add(model.trainingSet.id);
		if (model.hasValidationSets())
			for (Basket basket : model.getValidationSets())
				basketIdList.add(basket.id);

		Criteria c = Globals.session().createCriteria(Basket.class).add(Restrictions.in("id", basketIdList)).createAlias("entries", "entr")
				.createAlias("entr.ep", "exp").createAlias("exp.predicate", "pred").setProjection(projList);
		List<Object[]> predicatesUsed = c.list();

		for (Object[] row : predicatesUsed)
		{
			if ("-".equals(row[0].toString()))
				numIntervals += (Long) row[1];
			else if ("~".equals(row[0].toString()) || "~=".equals(row[0].toString()))
				numAppequals += (Long) row[1];
			else if (row[0].toString().contains(">") || row[0].toString().contains("<"))
				numGreaterless += (Long) row[1];
		}

		return new WebModel(getDefaultTemplate())
				.addParam("noChemaxonLicense","1")
				.addParam("numIntervals", numIntervals.toString()).addParam("numAppequals", numAppequals.toString())
				.addParam("numGreaterless", numGreaterless.toString()).setTemplate("modeller/configurators/data-preprocessing");
	}

	public void dataPreprocessingSubmit(HttpServletRequest request)
	{
		ModelProtocol protocol = model.attachment.getObject().protocol;

		DataPreprocessingParser.parseStandartizationUI(request, model.attachment.getObject().standartization, model.attachment.getObject().datahandling, protocol);
		
		this.currentPage = firstConfigurationStep;
	}

	/**
	 * Display the label weighting configuration dialog
	 * @return
	 */
	public WebModel labelWeighting()
	{
		LabelWeighting lw = getDefaultWeighting(model.trainingSet);
		return new WebModel(getDefaultTemplate()).addObject(lw).setTemplate("modeller/configurators/label-weighting");
	}

	/**
	 * Generate the default (equal) weighting schema based on a basket
	 * @param set
	 * @return
	 */
	@SuppressWarnings("unchecked")
	private LabelWeighting getDefaultWeighting(Basket set)
	{
		LabelWeighting weighting = new LabelWeighting();
		List<Object[]> rows = Globals
				.session()
				.createSQLQuery(
						"select Property.name as name, PropertyOption.name as opt from BasketEntry inner join ExperimentalProperty using (exp_property_id) inner join Property using(property_id) left join PropertyOption using (poption_id) where basket_id=:basketId group by Property.property_id, PropertyOption.poption_id")
				.setLong("basketId", set.id).list();
		for (Object[] row : rows)
			weighting.addClass((String) row[0], (String) row[1], 0.0);
		weighting.equalize();

		return weighting;
	}

	/**
	 * Configure the label weighting based on the user input
	 * @param request
	 */
	public void labelWeightingSubmit(HttpServletRequest request)
	{
		if (model.attachment.getObject().configuration instanceof CDSConfiguration)
			if (((CDSConfiguration) model.attachment.getObject().configuration).modelConfiguration != null &&
			((CDSConfiguration) model.attachment.getObject().configuration).modelConfiguration instanceof MultiLearningAbstractConfiguration
					)
			{
				MultiLearningAbstractConfiguration modelConf = (MultiLearningAbstractConfiguration)((CDSConfiguration) model.attachment.getObject().configuration).modelConfiguration;
				modelConf.labelWeighting = new LabelWeighting();
				if(request.getParameter("nomultilearning") != null)
					modelConf.noMultiLearning = true;

				model.implicitValues = request.getParameter("implicitValues");
				if(model.implicitValues != null && model.implicitValues.length() ==0)model.implicitValues = null;

				if(request.getParameter("global-weighting") != null)
					modelConf.labelWeighting.globalWeighting = true;
				if(request.getParameter("global-normalization") != null)
					modelConf.labelWeighting.globalNormalization = true;

				for (Object oName : request.getParameterMap().keySet())
				{
					String name = (String) oName;

					if (!name.contains("prop-"))
						continue;
					if (name.contains("$$$$"))
					{
						String[] parts = name.split("\\$\\$\\$\\$");
						modelConf.labelWeighting.addClass(resolvePropertyName(parts[0]).substring(5), parts[1].substring(6), new Double(request.getParameter(name)));
					}
					else
						modelConf.labelWeighting.addClass(resolvePropertyName(name.substring(5)), null, new Double(request.getParameter(name)));
				}

				for (Object oName : request.getParameterMap().keySet())
				{
					String name = (String) oName;

					if (name.startsWith("cm-"))
					{
						String[] parts = name.substring(3).split("--");
						if (request.getParameter("use-cm-" + parts[0]) != null)
						{
							modelConf.labelWeighting.getProperty(resolvePropertyName(parts[0])).setCostMatrix(parts[1], parts[2], new Double(request.getParameter(name)));
							modelConf.labelWeighting.useCostMatrix = true;
						}
					}
				}
				if(modelConf.labelWeighting.areWeightsIdentical()) {
					modelConf.labelWeighting.propertyWeights = null;
					if(!modelConf.labelWeighting.isGlobalWeighting())
						modelConf.labelWeighting = null;
				}
			}

		this.currentPage = "start";
	}

	//TODO This is a fix
	// Until correct Parent Property names will be shown 
	String resolvePropertyName(String name) {
		Property p = Repository.property.getProperty(name, false);
		if(p != null && p.parent != null) return p.parent.getName();
		return name;		
	}

}
