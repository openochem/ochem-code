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

import java.io.IOException;
import java.io.OutputStream;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.Hibernate;
import org.hibernate.HibernateException;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;
import org.hibernate.type.LongType;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.SessionVariable;
import qspr.ThreadScope;
import qspr.business.ModelOperation;
import qspr.business.PendingTaskPeer;
import qspr.business.Privileges;
import qspr.business.WebFilters;
import qspr.dao.Repository;
import qspr.entities.Alert;
import qspr.entities.Basket;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.ModelIdentity;
import qspr.entities.ModelMapping;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.ReadyModelAttachment;
import qspr.entities.User;
import qspr.export.ExportableColumn;
import qspr.export.ExportableModel;
import qspr.export.ExportableMolecule;
import qspr.export.ExportableSet;
import qspr.export.ExportableSetConfiguration;
import qspr.frontend.BrowserModel;
import qspr.frontend.ModelProfileData;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.metaserver.configurations.ADConfiguration;
import qspr.metaserver.configurations.LinearConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.PLSConfiguration;
import qspr.workflow.utils.QSPRConstants;
import qspr.modelling.ApplicabilityDomain;
import qspr.modelling.BestPredictionsPointSelector;
import qspr.modelling.ModelFactory;
import qspr.modelling.ModelProcessor;
import qspr.modelling.ModelStatistics;
import qspr.modelling.PointSelector;
import qspr.modelling.PointStatistics;
import qspr.modelling.ROCCurve;
import qspr.modelling.SetStatistics;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.applier.ModelApplierTaskProcessor;
import qspr.modelling.applier.PredictionResults;
import qspr.modelling.configurations.CDSModelData;
import qspr.util.AccessChecker;
import qspr.util.ExportThread;
import qspr.util.RWriter;
import qspr.util.RequestParser;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.useractions.EventFactory;
import com.eadmet.useractions.ModelCreationAction;
import com.eadmet.useractions.ModelCreationAction.ModelAction;
import com.eadmet.utils.NumericalValueStandardizer;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.mailer.Mailer;

@Controller
public class ModelController extends BrowserWrapper
{
	private static transient final Logger logger = LogManager.getLogger(ModelController.class);
	final static int POINTS_TO_SHOW  = 1000;

	public ModelController()
	{
		sessionRequired = true;
	}

	@NoLoginRequired
	@NoFrame
	public ModelAndView getFeaturedModelIDs(HttpServletRequest request, HttpServletResponse response) {
		@SuppressWarnings("unchecked")
		List<Long> ids = Globals.session().createSQLQuery("select published_id id from Model where featured_name is not null").addScalar("id", LongType.INSTANCE).list();
		return new ModelAndView("json", "object", ids);
	}

	public Model getRequestedModel()
	{
		if (assertParam("public_id"))
		{
			Model model = Repository.model.getByPublicId(getLongParam("public_id"));
			if (model != null && Globals.userSession() != null)
				Globals.userSession().visitedModelsIds.add(model.publicId);
			return model;
		}
		if (getLongParam("id") != null)
			return (Model) Globals.session().get(Model.class, getLongParam("id"));
		else
			return (Model) Globals.getSessionAttribute(SessionVariable.MODEL);
	}

	public ModelAndView getSizeSummary(HttpServletRequest request, HttpServletResponse response) throws Exception {
		Model m = getRequestedModel();
		Map<String, String> map = new HashMap<String, String>();	
		map.put("configuration", OCHEMUtils.getSizeBytes(m.attachment == null? 0: m.attachment.getDataLength()));
		map.put("model", OCHEMUtils.getSizeBytes(m.readyModelAttachment == null ? 0 : m.readyModelAttachment.getDataLength()));
		map.put("descriptors", OCHEMUtils.getSizeBytes(m.calcDescriptors == null ? 0: m.calcDescriptors.getDataLength()));
		map.put("statistics",OCHEMUtils.getSizeBytes(m.microattachment == null ? 0: m.microattachment.getDataLength()));

		return new ModelAndView("json", "object", map);
	}


	@NoFrame
	@NoLoginRequired
	public ModelAndView factSheet(HttpServletRequest request, HttpServletResponse response) {
		WebModel wm = new WebModel();
		wm.outerTemplate = "eadmet-outer";
		wm.templateName = "model/fact-sheet";

		Model model = getRequestedModel();
		wm.setObject(model);

		return wm.getModelAndView();
	}

	public ModelAndView markedForDeletionAction(HttpServletRequest request, HttpServletResponse response) throws HibernateException, Exception 
	{
		ModelOperation service = new ModelOperation();

		String action = getParam("action");
		if ("deleteMarkedModels".equals(action))
			service.processMarkedModels(true);
		else if ("keepMarkedModels".equals(action))
			service.processMarkedModels(false);
		else if ("deleteMarkedTasks".equals(action))
			PendingTaskPeer.processMarkedTasks(true);
		else if ("keepMarkedTasks".equals(action))
			PendingTaskPeer.processMarkedTasks(false);

		return new WebModel().getModelAndView();
	}

	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		WebFilters filters = formFilters(request);

		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);

		WebList list = new WebList();
		Criteria c = Globals.session().createCriteria(Model.class);

		// Filters
		if (assertParam("query"))
		{
			if (getParam("query").matches("[0-9]+"))
				c.add(Restrictions.or(Restrictions.like("name", "%" + getParam("query") + "%"), Restrictions.eq("publicId", getLongParam("query"))));
			else
				c.add(Restrictions.like("name", "%" + getParam("query") + "%"));
		}

		if (filters.has("publicId"))
			c.add(Restrictions.eq("publicId", getLongParam("publicId")));

		if (filters.has("awaiting-approval"))
		{
			filters.addFilterOverride("visibility", "public", null);
			c.add(Restrictions.eq("approved", false));
		}

		if (filters.has("status"))
		{
			if ("private".equals(filters.get("status")))
			{
				filters.addFilterOverride("visibility", "private", null);
			} else if ("awaiting".equals(filters.get("status")))
			{
				c.add(Restrictions.eq("approved", false));
				filters.addFilterOverride("visibility", "public", null);
			} else if ("approved".equals(filters.get("status")))
			{
				c.add(Restrictions.eq("approved", true));
				filters.addFilterOverride("visibility", "public", null);
			}
		}

		if (filters.has("featured"))
			c.add(Restrictions.isNotNull("featuredName")).add(Restrictions.eq("published", true));

		if (assertParam("toBeDeleted")) {
			c.add(Restrictions.isNotNull("deleteAfterDays")).add(Restrictions.eqOrIsNull("published", false));
		}

		if (assertParam("proquery"))
		{
			c.createAlias("modelMappings", "mm").createAlias("mm.property", "pro");
			c.add(Restrictions.like("pro.name", "%" + getParam("proquery") + "%"));
		}

		if (filters.has("property"))
		{
			c.createAlias("modelMappings", "mm");
			if (!filters.get("property").contains(","))
				c.add(Restrictions.eq("mm.property", Property.getById(filters.getLong("property"))));
			else
				c.createAlias("mm.property", "p").add(Restrictions.in("p.id", filters.getLongArray("property")));
		}

		if (assertParam("basket"))
		{
			Disjunction disj = Restrictions.disjunction();
			c.createAlias("validationSets", "vs", Criteria.LEFT_JOIN);
			disj.add(Restrictions.eq("trainingSet.id", getLongParam("basket")));
			disj.add(Restrictions.eq("validationSet.id", getLongParam("basket")));
			disj.add(Restrictions.eq("vs.id", getLongParam("basket")));
			c.add(disj);
		}

		if (assertParam("trainingSet"))
			c.add(Restrictions.eq("trainingSet.id", getLongParam("trainingSet")));

		if (assertParam("article"))
		{
			if (getParam("article").matches("a.*|A.*"))
				c.add(Restrictions.eq("article.id", Long.valueOf(getParam("article").substring(1))));
			else
				c.add(Restrictions.eq("article.id", getLongParam("article")));
		}

		String visibility = filters.get("visibility");
		if(visibility == null)visibility="null";
		switch(visibility) {
		case "publishedAll" :
		case "public" : c.add(Restrictions.eq("published", Boolean.TRUE)); break;
		case "private" : c.add(Restrictions.eq("published", Boolean.FALSE)); break;
		default: break;  // no restriction for other categories
		}

		if(!Globals.isOCHEMDeveloper())
			visibility = "all";

		switch(visibility) {
		case "all" : // these three only within the group
		case "public" :
		case "private" : 
			Model.addAuthRestrictions(c, assertParam("group"), filters.has("awaiting-approval") || "awaiting".equals(getParam("status"))); break; // only for this user
		default: break; // otherwise no restriction
		}

		// only models and not pending tasks !
		c.add(Restrictions.isNull("taskId"));

		if (assertParam("introducer"))
			c.add(Restrictions.eq("sess.user", User.getByString(getParam("introducer"))));

		if (assertParam("template"))
		{
			c.createAlias("template", "t");
			c.add(Restrictions.eq("t.name", getParam("template")));
		}

		if (assertParam("order"))
		{
			if (getParam("order").equals("access"))
				c.addOrder(Order.desc("lastAccess"));
			else if (getParam("order").equals("modification"))
				c.addOrder(Order.desc("lastModification"));
			else
				c.addOrder(Order.desc("id")); //Creation time
		}

		list.useEntity(Model.class).loadDistinctFromCriteria(c, getPageNum(), getPageSize(15));

		// Add checkboxes
		if (applier != null)
		{
			for (Object elem : list.list)
			{
				Model m = (Model) elem;
				m.selected = (applier.indexOf(m) != null);
			}
		}
		return new BrowserModel().setFilters(filters).setObject(list).getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		String action = request.getParameter("action");
		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);

		Long id = Long.valueOf(request.getParameter("id"));
		Model model = (Model) Globals.session().get(Model.class, id);
		Privileges privileges = model.getPrivileges(request);

		if (!privileges.canView)
			throw new UserFriendlyException("Not sufficient privileges");

		if (action.equals("publish"))
		{
			if (!privileges.canEdit)
				throw new UserFriendlyException("Not sufficient privileges");

			long articleId = Long.valueOf(request.getParameter("article"));

			model = (Model) Globals.session().get(Model.class, id); // re-fetching the model
			model.publish(articleId);
			Globals.session().saveOrUpdate(model);

			model = (Model) Globals.session().get(Model.class, id); // re-fetching the model
			model.publishModelPredictions(articleId); 
			Mailer.notifyAdmins("A published model is awaiting approval", "A model " + model.name + " with public ID " +model.publicId + " has been published by user " + model.session.user.login + ". It is not yet publicly listed - it is awaiting for you approval.\nUse your favourite SQL query tool to approve the model.");

			EventFactory.document("Model publish", null, "has published a model for " + model.modelMappings.get(0).property.getName());

			return redirect("model/profile.do?id=" + model.id);
		}
		else if (action.equals("delete"))
		{
			if (!privileges.canEdit)
				throw new UserFriendlyException("Not sufficient privileges");
			model.delete();
			EventFactory.document("Model deletion", new ModelCreationAction(model, ModelAction.DELETE), null);
		}
		else if ("togglebasket".equals(action))
		{
			if (applier != null)
			{
				model.markAccessed();
				Globals.session().saveOrUpdate(model);
				Integer index = applier.indexOf(model);
				if (index != null)
					applier.removeModel(model);
				else
					applier.addModel(model);
				model.selected = new Boolean(index == null);
			}
		}
		else if ("rename".equals(action))
		{
			model.name = getParam("name");
			Globals.session().saveOrUpdate(model);
			logger.info("Renaming...");
		}

		return new WebModel().getModelAndView();
	}

	public ModelAndView select(HttpServletRequest request, HttpServletResponse response)
	{
		ModelApplier applier = new ModelApplier();
		Globals.setSessionAttribute(SessionVariable.MODEL_APPLIER, applier);
		String template = request.getParameter("browser") != null ? "model/browser" : "model/select";
		return new WebModel().setTemplate(template).getModelAndView();
	}

	// Create a clone of a model
	public ModelAndView createModelCopy(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Model model = getRequestedModel();

		Privileges privileges = model.getPrivileges(request);
		if (!privileges.canEdit)
			throw new UserFriendlyException("You are not authorized to copy this model");

		Model copy = new Model();
		copy.session = Globals.userSession();
		if (model.calcDescriptors != null)
			copy.calcDescriptors = model.calcDescriptors.getCopy();
		copy.configurationXml = model.configurationXml;
		copy.attachment = model.attachment.getCopy();
		if(model.readyModelAttachment != null)copy.readyModelAttachment = model.readyModelAttachment.getCopy();
		copy.defaultTaskPriority = model.defaultTaskPriority;
		copy.fullEquation = model.fullEquation;
		copy.isStatisticsCalculated = model.isStatisticsCalculated;
		copy.lastAccess = copy.dateCreated = copy.lastModification = new Timestamp(Calendar.getInstance().getTimeInMillis());
		copy.description = model.description;
		copy.microattachment = model.microattachment.getCopy();
		copy.name = "Copy of " + model.name;
		copy.predictionScenario = model.predictionScenario;
		copy.trainingSet = model.trainingSet;
		for (Basket b : model.getValidationSets())
			copy.addValidationSet(b);
		copy.template = model.template;
		copy.timeToComplete = model.timeToComplete;

		Globals.session().save(copy);

		for (ModelMapping mm : model.modelMappings)
		{
			ModelMapping mmCopy = new ModelMapping();
			mmCopy.model = copy;
			mmCopy.property = mm.property;
			mmCopy.statisticsOriginal = mm.statisticsOriginal.getCopy();
			mmCopy.statisticsRecalculated = mm.statisticsRecalculated.getCopy();
			mmCopy.unit = mm.unit;
			mmCopy._class = mm._class;
			Globals.session().save(mmCopy);
		}

		return redirect("model/profile.do?id=" + copy.id);
	}

	public ModelAndView iterationsPlot(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Thread.currentThread().setContextClassLoader(ModelController.class.getClassLoader());
		Model model = getRequestedModel();

		Privileges privileges = model.getPrivileges(request);
		if (!privileges.canView)
			throw new UserFriendlyException("Not sufficient privileges");

		List<ROCCurve> plot = new ArrayList<ROCCurve>();

		if (model.getModelData(true) instanceof CDSModelData)
		{
			CDSModelData cdsData = (CDSModelData) model.getModelData(true);
			if(((ModelAbstractConfiguration) cdsData.methodSpecificData).iterations != null) {
				DataTable iterations = ((ModelAbstractConfiguration) cdsData.methodSpecificData).iterations;
				if(iterations.getRowsSize() > 3) {
					iterations.deleteRow(0);  // first and last rows can be garbage...
					iterations.deleteRow(iterations.getRowsSize()-1);
					int valid = iterations.getColumnsSize()-1; // last column contains validation set performance
					for(int i=0;i<=valid;i++) 
						plot.add(new ROCCurve());
					plot.get(0).setId = "Training";
					plot.get(1).setId = "Validation";
					plot.get(valid).setId = "Early stop";

					if(iterations != null && iterations.getRowsSize() > 0) {
						int rowMinTest = iterations.getMinValueRowId(2, iterations.getRowsSize()/100 > 1?iterations.getRowsSize()/100:1);
						double minValValue = (double)iterations.getValue(rowMinTest,2);
						for(int i=0;i<iterations.getRowsSize();i++) {
							double iter = (double)iterations.getValue(i, 0);
							for(int j=0;j<valid;j++) {
								Double v = (Double)iterations.getValue(i, j + 1);
								if(v == null)continue;
								plot.get(j).addPoint(iter, v, (long)i);
							}
						}
						for(int i = 0;i <= rowMinTest; i++) // add the early stopping curve
							plot.get(valid).addPoint((Double)iterations.getValue(i, 0), minValValue, (long)i);
					}
				}
			}
		}

		return new WebModel().setList(plot).getModelAndView();
	}

	public ModelAndView rocCurve(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Thread.currentThread().setContextClassLoader(ModelController.class.getClassLoader());
		Model model = getRequestedModel();
		ModelMapping selectedMapping = null;
		if (model == null || assertParam("mapping_id"))
		{
			if (!assertParam("mapping_id"))
				throw new UserFriendlyException("The model you a trying to access is not available or has been deleted");
			selectedMapping = (ModelMapping) Globals.session().get(ModelMapping.class, getLongParam("mapping_id"));
			model = selectedMapping.model;
		}
		if (selectedMapping == null)
			selectedMapping = model.modelMappings.get(0);

		Privileges privileges = model.getPrivileges(request);
		if (!privileges.canView)
			throw new UserFriendlyException("Not sufficient privileges");

		ModelStatistics ms = (ModelStatistics) selectedMapping.statisticsRecalculated.getObject();
		List<ROCCurve> res = new ArrayList<ROCCurve>();
		Long validationSetId = getLongParam(QSPRConstants.VALIDATION);
		for (int set = 0; set < ms.sets.size(); set++)
		{
			SetStatistics ss = ms.sets.get(set);
			ROCCurve curve = ss.getROCCurve();
			if (ss.setId.startsWith(QSPRConstants.VALIDATION))
			{
				if (ss.basketId == null || !ss.basketId.equals(validationSetId))
					continue;
				curve.setId = QSPRConstants.VALIDATION;
			}
			res.add(curve);
		}
		return new WebModel().setList(res).getModelAndView();
	}

	public ModelAndView profile(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ThreadScope.get().considerPredicatesInStatistics = assertParam("consider-predicates");
		Globals.setMarshallingOption(MarshallingOption.MODEL_SET_SIZE);
		Thread.currentThread().setContextClassLoader(ModelController.class.getClassLoader());
		List<Object> list = new ArrayList<Object>();


		Model model = getRequestedModel();
		ModelMapping selectedMapping = null;
		String validationSetId = getParam(QSPRConstants.VALIDATION);

		if (model == null || assertParam("mapping_id"))
		{
			if (!assertParam("mapping_id"))
				throw new UserFriendlyException("The model you a trying to access is not available or has been deleted");
			selectedMapping = (ModelMapping) Globals.session().get(ModelMapping.class, getLongParam("mapping_id"));
			model = selectedMapping.model;
		}

		List<ModelMapping> modelMList = model.modelMappings;
		if (modelMList.size() == 1)
			selectedMapping = model.modelMappings.get(0);

		Privileges privileges = model.getPrivileges(request);
		if (!privileges.canView)
			throw new UserFriendlyException("Not sufficient privileges");

		if (privileges.canEdit)
		{

			if (assertParam("delete-validation-set")) // for debugging only
			{
				Long setId = "all".equals(getParam("delete-validation-set")) ? null : getLongParam("delete-validation-set");
				model.deleteValidationSet(setId);
			}

			if (assertParam("add-validation-set"))
			{
				ModelApplier applier = getModelApplier();
				for (ModelApplierTaskProcessor mTask : applier.modelTasks)
					if (mTask.model.id.equals(model.id))
					{
						mTask.addValidationSet();
					}

				return redirect("model/profile.do?id=" + model.id + "&render-mode=popup");
			}

			if (assertParam("restore-basket"))
				if (model.restoreModifiedBasketEntries("recalculated".equals(getParam("restore-basket"))))
					return redirect("model/profile.do?id=" + model.id + "&mapping_id=" + selectedMapping.id + "&render-mode=popup");

			if (assertParam("exclude-duplicates"))
				if (model.excludeDuplicates())
					return redirect("model/profile.do?id=" + model.id + "&mapping_id=" + 
							(selectedMapping == null? model.modelMappings.get(0).id:selectedMapping.id)
							+ "&render-mode=popup");

			if (assertParam("storeRecalculatedModel"))
			{
				// Replace the original model with the recalculated one
				model.replaceOriginalWithRecalculated();
			}
			else if (assertParam("deleteRecalculatedModel"))
			{
				// Delete the recalculated model
				model.setModelData(model.getModelData(false), true);
				model.attachment.updateObject();
				model.readyModelAttachment.updateObject();
				for (ModelMapping mm : model.modelMappings)
				{
					ModelStatistics msOriginal = (ModelStatistics) mm.statisticsOriginal.getObject();
					mm.statisticsRecalculated.setObject(msOriginal);
					mm.statisticsOriginal.updateObject();
					Globals.session().saveOrUpdate(mm);
				}
				if (model.restoreModifiedBasketEntries(false))
					return redirect("model/profile.do?id=" + model.id + "&mapping_id=" + selectedMapping.id + "&render-mode=popup");

				Globals.session().saveOrUpdate(model);
			}
		}

		if (selectedMapping != null)
		{
			Hibernate.initialize(selectedMapping.statisticsOriginal);
			Hibernate.initialize(selectedMapping.statisticsRecalculated);

			if (selectedMapping.statisticsOriginal == null)
			{
				ModelProcessor processor = ModelFactory.getProcessor(model.template);
				processor.model = model;
				processor.wrapped();
				Globals.session().saveOrUpdate(model);
			}

			logger.info("Recalculating statistics");
			ModelStatistics msFirst = null, msSecond = null;

			PointSelector pointSelector = getPointSelector(request);

			if (!selectedMapping.statisticsOriginal.getAttachmentReference().equals(selectedMapping.statisticsRecalculated.getAttachmentReference()))
			{
				msFirst = (ModelStatistics) selectedMapping.statisticsOriginal.getObject();
				list.add(msFirst);
				if (((ModelStatistics) selectedMapping.statisticsRecalculated.getObject()).containsPredictions)
				{
					msSecond = (ModelStatistics) selectedMapping.statisticsRecalculated
							.getObject();
					msSecond.validationSetId = validationSetId;
					msSecond.actualizeStatistics(selectedMapping);
					msSecond.recalculateStatistics(selectedMapping, pointSelector);
					list.add(msSecond);
				}
			}
			else
			{
				msFirst = (ModelStatistics) selectedMapping.statisticsRecalculated.getObject();
				list.add(msFirst);
			}

			msFirst.validationSetId = validationSetId;
			msFirst.actualizeStatistics(selectedMapping);
			msFirst.recalculateStatistics(selectedMapping, pointSelector);

			ModelStatistics modelStatistics = msSecond != null ? msSecond : msFirst;

			Pattern p = Pattern.compile("[\\u0009\\u000A\\u000D\u0020-\\uD7FF\\uE000-\\uFFFD\\u10000-\\u10FFF]+"); 
			for (int i=0; i<modelStatistics.sets.get(0).points.size(); i++)
				if (modelStatistics.sets.get(0).points.get(i).error != null)
					if (!p.matcher(modelStatistics.sets.get(0).points.get(i).error).matches())
						modelStatistics.sets.get(0).points.get(i).error = "An error was found while processing this molecule. In addition, the error message contained non-UTF characters and was removed.";

			if (modelStatistics.sets.get(0).distancesToModel.size() > 0)
			{

				int dmNum = assertParam("dm") ? getIntParam("dm") : 0;
				String dmName = assertParam("dmname") ? getParam("dmname")
						: modelStatistics.sets.get(0).distancesToModel.get(dmNum);

				logger.info("Calculating applicability domain: " + dmName);

				if (dmName != null)
				{
					for (String dm : modelStatistics.sets.get(0).distancesToModel)
					{
						if (dm.equals(dmName) || dmName.equals("all"))
						{
							ApplicabilityDomain ad = new ApplicabilityDomain();
							ad.showNegativeValues = assertParam("show-negative");
							//ad.useStandardizedResiduals = dmName.equals("Leverage") && AccessChecker.isCadaster(Globals.userSession().user);
							ad.showNegativeValues |= ad.useStandardizedResiduals;

							ADConfiguration adConfiguration = new ADConfiguration();
							if (assertParam("wsize"))
								adConfiguration.windowSizeInPercent = Double.valueOf(getParam("wsize")) / 100;
							if (assertParam("averaging"))
								adConfiguration.averagingType = getParam("averaging");
							ad.setModel(selectedMapping, dm, adConfiguration);
							// Keep necessary data for X-axis
							for (SetStatistics ss : modelStatistics.sets)
							{
								if ("percents".equals(getParam("xaxis")))
								{
									if (ss.adConfiguration != null)
										ss.adConfiguration.intervals = null;
								}
								else
								{
									if (ss.adConfiguration != null)
										ss.adConfiguration.percents = null;
								}
							}
							list.add(ad);
						}
					}
				}
				logger.info("Calculating AD finished");
			}

			if (selectedMapping.property.isQualitative())
				list.addAll(fillInOptions(selectedMapping));
		}

		if (model.id == null && assertParam("save"))
		{
			model.updateDescription();
			Globals.session().save(model);
			for (ModelMapping mm : model.modelMappings)
				Globals.session().saveOrUpdate(mm);
		}

		if (assertParam("updateDescription"))
		{
			model.updateDescription();
			Globals.session().saveOrUpdate(model);
		}

		if (model.size == null || model.size <= 300 || assertParam("updateDescription"))
			model.recalculateModelSize();

		// should the model profile ever be displayed,
		// remove the model immediately from pending list
		if (!assertParam("save"))
		{
			model.markAccessed();
			Globals.session().saveOrUpdate(model);
			Globals.session().flush();
		}

		if (selectedMapping != null)
		{
			model.modelMappings.clear();
			model.modelMappings.add(selectedMapping);
			Globals.session().evict(model);
		}

		return new WebModel(model).setList(list).setTemplate(
				selectedMapping != null ? "model/single-profile" : "model/multi-profile").setRenderMode("popup")
				.getModelAndView();
	}


	public ModelAndView exportDescriptors(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		final ExportableSet eData = new ExportableSet();
		eData.removeColumn(ExportableColumn.APPLICABILITY_DOMAIN); 
		eData.removeColumn(ExportableColumn.PREDICTED_VALUE); 
		eData.removeColumn(ExportableColumn.EXP_VALUE);
		eData.removeColumn(ExportableColumn.EXP_VALUE_CONVERTED); 
		eData.removeColumn(ExportableColumn.DM_VALUE);
		eData.removeColumn(ExportableColumn.ACCURACY);

		if (assertParam("submit"))
		{
			final String format = getParam("format");
			final long modelId = getLongParam("model");

			ExportThread eThread = new ExportThread(format, ExportableSetConfiguration.configureFromDialog(request))
			{
				@Override
				public void generateData() throws Exception
				{
					DataTable dtDescriptors;
					Model model;
					model = (Model) Globals.session().get(Model.class, modelId);
					eData.addModel(model);
					eData.setDescriptors(dtDescriptors = model.getCalculatedDescriptors());

					dtDescriptors.reset();
					System.out.println(dtDescriptors.toStringColumns());
					while (dtDescriptors.nextRow())
					{
						if (eData.exportableMolecules.size() % 1000 == 0 && eData.exportableMolecules.size() > 0)
							Globals.restartAllTransactions(true);
						if (eData.exportableMolecules.size() % 100 == 0)
							setStatus("Preparing item "+eData.exportableMolecules.size());
						ExportableMolecule eMol = new ExportableMolecule();
						eData.addMolecule(eMol);
						ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, (Long)dtDescriptors.getCurrentRow().getAttachment(QSPRConstants.RECORD_ID_ATTACHMENT));
						eMol.setExperimentalProperty(ep);
						eMol.descriptors = dtDescriptors.getCurrentRow();
					}

					eData.supplementaryData.put("Configuration", model.configurationXml);
					setFileName("descriptors_"+model.name);
				}
			};
			eThread.start();
			Thread.sleep(100);
			return redirect("longoperations/operationWaitingScreen.do?operation-id="+eThread.operationID);
		}
		else
		{
			return new WebModel(eData).setTemplate("export").getModelAndView();
		}
	}

	public ModelAndView setoverview(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		// Retrieve all mappings, not models
		Criteria criteria = Globals.session().createCriteria(ModelMapping.class);

		Criteria c = criteria.createCriteria("model");
		c.add(Restrictions.or(Restrictions.eq("trainingSet.id", getLongParam("basket")), Restrictions.eq(
				"validationSet.id", getLongParam("basket"))));
		Model.addAuthRestrictions(c, assertParam("group"), false);
		c.add(Restrictions.isNull("taskId"));

		return new WebModel().setList(criteria.list()).setTemplate("model/basket-models").getModelAndView();
	}

	public ModelAndView loadAd(HttpServletRequest request, HttpServletResponse response)
	{
		Model model = getRequestedModel();
		ModelMapping selectedMapping = null;

		if (model == null || assertParam("mapping_id"))
		{
			if (!assertParam("mapping_id"))
				throw new UserFriendlyException("The model you a trying to access is not available or has been deleted");
			selectedMapping = (ModelMapping) Globals.session().get(ModelMapping.class, getLongParam("mapping_id"));
			model = selectedMapping.model;
		}

		Privileges privileges = model.getPrivileges(request);
		if (!privileges.canView)
			throw new UserFriendlyException("Not sufficient privileges");

		ModelStatistics ms = (ModelStatistics) selectedMapping.statisticsRecalculated.getObject();
		SetStatistics ss = ms.sets.get(0);

		ADConfiguration adConf = new ADConfiguration();
		adConf.averagingType = "cumulative";

		String dmName = assertParam("dmname") ? getParam("dmname") : ms.sets.get(0).distancesToModel.get(0);

		adConf = ss.getADConfiguration(dmName, null, selectedMapping, adConf);

		return new WebModel(adConf).getModelAndView();
	}

	public ModelAndView loadstatistics(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Globals.setMarshallingOption(MarshallingOption.MODEL_SET_SIZE);

		if (assertParam("articles"))
			Globals.setMarshallingOption(MarshallingOption.ARTICLE_IN_MODEL_DOT);

		List<Object> others = new ArrayList<Object>();

		ModelMapping mm = null;

		if (assertParam("id"))
		{
			Model model = (Model) Globals.session().get(Model.class, getLongParam("id"));
			mm = model.modelMappings.get(0);
			others.add(model);
		}

		if (assertParam("mm_id"))
		{
			mm = (ModelMapping) Globals.session().get(ModelMapping.class, getLongParam("mm_id"));
			others.add(mm);
		}

		ModelStatistics stats = (ModelStatistics) mm.statisticsRecalculated.getObject();
		//stats.actualizeStatistics(mm);
		//stats.recalculateStatistics(mm);

		if (mm.property.isQualitative())
		{
			for (SetStatistics set : stats.sets)
				set.points = null;
			others.addAll(fillInOptions(mm));
		}else{
			for (SetStatistics set : stats.sets)
				if(set.points.size() > POINTS_TO_SHOW) {
					Collections.shuffle(set.points);
					set.points = new ArrayList<PointStatistics>(set.points.subList(0, POINTS_TO_SHOW));
				}
		}

		mm.model.markAccessed();
		Globals.session().saveOrUpdate(mm.model);

		return new WebModel(stats).setList(others).getModelAndView();
	}

	/**
	 * Filling confusion tables or added empty values if something does not exits
	 * @param selectedMapping
	 * @return
	 */

	List<PropertyOption> fillInOptions(ModelMapping selectedMapping) {
		List<PropertyOption> list = new ArrayList<PropertyOption>();

		if (selectedMapping.property.isQualitative())
		{
			List<PropertyOption> options = new ArrayList<PropertyOption>();
			for (Long optionId : selectedMapping.model.attachment.getObject().optionsMapping.keySet())
			{
				PropertyOption po = (PropertyOption) Globals.session().get(PropertyOption.class, optionId);
				if (po.property.equals(selectedMapping.property))
				{
					int _class = selectedMapping.model.getMappedOption(optionId).intValue();
					while (options.size() <= _class)
						options.add(new PropertyOption());
					options.set(_class, po);
				}
			}

			list.addAll(options);
		}
		return list;
	}	

	public ModelAndView exportModelAsR(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Model model = (Model) Globals.session().get(Model.class, getLongParam("id"));

		model.trainingSet.getPrivileges().requestView();

		ModelStatistics stats = (ModelStatistics) model.modelMappings.get(0).statisticsRecalculated.getObject();

		response.setContentType("application/R");
		response.setHeader("Content-Disposition", "attachment; filename=\"model_" + model.name.replaceAll("[ \\.]", "_")
		+ ".R\"");

		RWriter writer = new RWriter(response.getOutputStream());
		writer.variablePrefix = model.name.replaceAll("\\s", ".") + "$";

		for (int set = 0; set < stats.sets.size(); set++)
			for (int pnum = 0; pnum < stats.sets.get(set).points.size(); pnum++)
			{
				PointStatistics point = stats.sets.get(set).points.get(pnum);
				if (pnum % 100 == 0)
					logger.info("" + pnum + " records processed");
				if (point.error != null)
					continue;
				ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class,
						point.id);
				// logger.info(point.real);
				// logger.info(point.distancesToModel);
				writer.addValue("set" + set + "$dm", "" + stats.sets.get(set).getDM(pnum, null));
				// writer.addValue("set" + set + "$num", ep.artMolId);
				writer.addValue("set" + set + "$epid", "" + point.id);
				writer.addValue("set" + set + "$mw", "" + ep.molecule.molWeight);
				writer.addValue("set" + set + "$mp1", "" + ep.molecule.mapping1.id);
				writer.addValue("set" + set + "$mp2", "" + ep.molecule.mapping2.id);
				writer.addValue("set" + set + "$obs.value", "" + point.real);
				writer.addValue("set" + set + "$pred.value", "" + point.predicted);
			}

		DataTable dtDescriptors = model.getCalculatedDescriptors();
		dtDescriptors.reset();
		while (dtDescriptors.nextRow())
			writer.addValue("descriptors", "" + dtDescriptors.getValue());

		writer.writeAndClose();

		return null;
	}

	public ModelAndView browseRecords(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ModelMapping mm = (ModelMapping) Globals.session().get(ModelMapping.class, getLongParam("mm_id"));
		String setId = getParam("set");
		ModelStatistics ms = (ModelStatistics) (assertParam("recalculated") ? mm.statisticsRecalculated.getObject() : mm.statisticsOriginal.getObject());
		if (assertParam("validationSetNum"))
			ms.validationSetId = getParam("validationSetNum");
		ms.actualizeStatistics(mm);
		SetStatistics ss = ms.getSetStatistics(setId);

		if(ss==null && setId.equals("validation"))ss = ms.getSetStatistics(setId + 0); // bug since 0 set has id 0

		// A very unefficient and not nice way. May completely fail for large datasets. Think of optimisation/refactoring later / Midnighter on Jun 16, 2011
		// An old problem - model records' IDs are stored inside of a serializable object and therefore it is tricky to query them. See also Bug 1809
		if (ss==null || ss.points.size() < 1000)
		{
			StringBuffer buff = new StringBuffer();
			if(ss != null)for (PointStatistics ps : ss.points)
			{
				buff.append(ps.id);
				buff.append(",");
			}
			buff.append(0);

			return redirect("epbrowser/show.do?id=" + buff.toString());
		}
		else
		{
			// Temporary heuristics. Should be implemented in a more robust way later / Midnighter on Jun 16, 2011
			if (ss.basketId != null && !QSPRConstants.EXCLUDED.equals(setId))
				return redirect("epbrowser/show.do?basket-select=" + ss.basketId + "&excluded=-" + mm.model.id + "&property=" + mm.property.id);
			else if (QSPRConstants.TRAINING.equals(setId))
				return redirect("epbrowser/show.do?basket-select=" + mm.model.trainingSet.id + "&excluded=-" + mm.model.id + "&property=" + mm.property.id);
			else if (QSPRConstants.VALIDATION.equals(setId))
				return redirect("epbrowser/show.do?basket-select=" + mm.model.getValidationSets().get(0).id + "&excluded=-" + mm.model.id + "&property=" + mm.property.id);
			else if (QSPRConstants.EXCLUDED.equals(setId))
				return redirect("epbrowser/show.do?basket-select=" + mm.model.trainingSet.id + "&excluded=" + mm.model.id + "&property=" + mm.property.id);
			else
				throw new UserFriendlyException("Sorry, the results cannot be displayed for technical reasons");
		}
	}

	public ModelAndView exportModel(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		final ExportableSet eData = new ExportableSet();
		eData.removeColumn(ExportableColumn.ACCURACY);
		eData.uncheckColumn(ExportableColumn.DESCRIPTORS);
		final Model model = (Model) Globals.session().get(Model.class, getLongParam("id"));

		model.trainingSet.getPrivileges().requestView();

		for (ModelMapping mm : model.modelMappings)
			eData.properties.add(mm.property);
		// If a model contains restricted access descriptors, disallow to export descriptors at all
		if (!AccessChecker.canExportDescriptors(Globals.userSession().user, model.configurationXml == null ? model.description : model.configurationXml.toLowerCase()))
			eData.removeColumn(ExportableColumn.DESCRIPTORS);
		else
			eData.removeColumn(ExportableColumn.DESCRIPTORSNAMES);

		if (assertParam("submit"))
		{
			eData.configure(ExportableSetConfiguration.configureFromDialog(request));
			eData.supplementaryData.put("Model configuration", model.configurationXml);

			final long modelId = Long.valueOf(getLongParam("id"));
			final String format = getParam("format");

			ExportThread wp = new ExportThread(format, ExportableSetConfiguration.configureFromDialog(request))
			{
				@Override
				public void generateData() throws Exception
				{
					Model model = (Model) Globals.session().get(Model.class, modelId);
					Hibernate.initialize(model.attachment);
					ExportableMolecule eMol;
					DataTable dtDescriptors = null;

					int sheetNum = 0;
					int totalSize = 0;
					for (ModelMapping mm : model.modelMappings)
						totalSize += ((ModelStatistics) mm.statisticsRecalculated.getObject()).getRowsSize();

					for (ModelMapping mm : model.modelMappings)
					{
						mm = (ModelMapping) Globals.session().get(ModelMapping.class, mm.id);
						// create separate worksheet for each property
						ModelStatistics statistics = (ModelStatistics) mm.statisticsRecalculated.getObject();

						// FIXME: Add authorization block to the universal export
						if (model.calcDescriptors != null && 
								(eData.selectedColumns.contains(ExportableColumn.DESCRIPTORS) ||
										eData.selectedColumns.contains(ExportableColumn.DESCRIPTORSNAMES))
								)
						{
							eData.setDescriptors(dtDescriptors = model.getCalculatedDescriptors());
							dtDescriptors.reset();
						}

						for (int setNum = 0; setNum < statistics.sets.size(); setNum++)
						{
							SetStatistics setStatistics = statistics.sets.get(setNum);
							setStatus("Preparing dataset " + setStatistics.setId);

							if (setStatistics.points.size() == 0)
								continue;

							String propertyName = mm.property.getName();
							eData.useSheet("" + (++sheetNum) + "_" + propertyName + "_" + setStatistics.setId);

							long lastFetchedEntry = -1;
							for (int entryNum = 0; entryNum < setStatistics.points.size(); entryNum++)
							{
								if (lastFetchedEntry < entryNum)
								{
									Globals.restartAllTransactions(true);
									List<Long> epIDs = new ArrayList<Long>();
									for (int i = entryNum; i < setStatistics.points.size() && i <= entryNum + 200; i++)
										epIDs.add(setStatistics.points.get(i).id);
									//setStatus("Preloaded " + entryNum + " entries out of " + setStatistics.points.size() + " in set " + setStatistics.setId);
									Globals.session().createCriteria(ExperimentalProperty.class)
									.createAlias("molecule", "mol")
									.createAlias("mol.mapping1", "m1")
									.createAlias("mol.mapping2", "m2")
									.add(Restrictions.in("id", epIDs)).list();
									lastFetchedEntry = entryNum + 200;
								}

								eData.addMolecule(eMol = new ExportableMolecule());
								PointStatistics point = setStatistics.points.get(entryNum);
								Long id = point.id;
								AbstractDataRow descriptorsRow = null;
								// Find the appropriate row in the descriptors table
								if (dtDescriptors != null)
								{
									dtDescriptors.nextRow();
									try
									{
										while (!id.equals(dtDescriptors.getCurrentRow().getAttachment(QSPRConstants.RECORD_ID_ATTACHMENT)))
											dtDescriptors.forceNextRow();
										descriptorsRow = dtDescriptors.getCurrentRow();
									}
									catch (Exception e)
									{

									}
								}

								ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, id);

								if (ep == null) {
									eMol.setExperimentalProperty(point, mm.property, mm.unit, model.attachment.getObject());
									if(point.virtual == null)
										eMol.error = "Record " + id + " has been deleted from the database";
								}
								else
									eMol.setExperimentalProperty(ep);

								if( (ep == null || ep.isDeleted()) && !Globals.isOCHEMDeveloper())  // deleting information for deleted (hidden) molecules
									eMol.clearSensitive();


								if (eMol.error == null && point.error != null)
									eMol.error = point.error;

								if(point.error == null)
								{
									double predicted = point.predicted;
									if (mm.property.isQualitative())
										eMol.setPrediction(mm, model.attachment.getObject().getOptionFromPrediction(predicted, mm.property).name, predicted);
									else
										eMol.setPrediction(mm, predicted);
									// Distances to model
									int dmNum = 0;
									for (String dm : setStatistics.distancesToModel)
									{
										eMol.setDM(mm, dm, point.distancesToModel.get(dmNum));
										dmNum++;
									}

									// Descriptors
									if (dtDescriptors != null)
									{
										eMol.descriptors = descriptorsRow;
									}
								}

								if (eData.exportableMolecules.size() % 100 == 0)
									setStatus("Preparing item " + eData.exportableMolecules.size() + " out of " + totalSize);
							}
						}
					}

					ModelAbstractConfiguration modelConfiguration = null; 
					model = (Model)Globals.session().get(Model.class, model.id); //Re-fetch model, otherwise some fields of it get detached

					if (model.getModelData(true) instanceof CDSModelData)
					{
						CDSModelData cdsData = (CDSModelData) model.getModelData(true);
						modelConfiguration = (ModelAbstractConfiguration) cdsData.methodSpecificData;
						if (modelConfiguration == null || modelConfiguration.iterations == null || modelConfiguration.iterations.getRowsSize() == 0)
							modelConfiguration.iterations = new DataTable(true);
						String name = model.name + " PublicId: " + model.publicId;
						modelConfiguration.iterations.addColumn(name);
						int iter = modelConfiguration.iterations.getMinValueRowId(2,1);
						if(iter > 0) {
							modelConfiguration.iterations.setValue(0, name, (double)iter);
							modelConfiguration.iterations.setValue(1, name, modelConfiguration.iterations.getValue(iter, 0));
							if(modelConfiguration.iterations.getValue(iter, 1) != null)
								modelConfiguration.iterations.setValue(2, name, modelConfiguration.iterations.getValue(iter, 1));
							if(modelConfiguration.iterations.getValue(iter, 2) != null)
								modelConfiguration.iterations.setValue(3, name, modelConfiguration.iterations.getValue(iter, 2));
						}
						if(modelConfiguration.versionOCHEM != null)
							modelConfiguration.iterations.addColumn(modelConfiguration.versionOCHEM);
						eData.supplementaryData.put("Info", modelConfiguration.iterations);
					}

					// Export of a linear model equation
					if (modelConfiguration instanceof LinearConfiguration ) 
					{	
						LinearConfiguration linearData = (LinearConfiguration) modelConfiguration;
						if(model.calcDescriptors != null) {
							eData.supplementaryData.put("Equation", 
									linearData.getEquation(model.getCalculatedDescriptors().getColumns(), false));
							eData.supplementaryData.put("Normalised equation", 
									linearData.getEquation(model.getCalculatedDescriptors().getColumns(), true));

							if (modelConfiguration instanceof PLSConfiguration)
							{
								PLSConfiguration plsData = (PLSConfiguration) modelConfiguration;
								if (plsData.lvMatrix != null)
								{
									Map<String, String> lvEquations = new LinkedHashMap<String, String>();
									lvEquations.put("Bias", plsData.coefficient.get(0).toString());
									for (int i = 0; i < plsData.lvMatrix.length; i++)
									{
										StringBuilder eq = new StringBuilder();
										for (int k = 0; k < plsData.lvMatrix[i].length; k++)
										{
											String val = NumericalValueStandardizer.getSignificantDigitsStr(Math.abs(plsData.lvMatrix[i][k]), 3);
											eq.append(plsData.lvMatrix[i][k] >= 0 ? " + " : " - ");
											eq.append(val);
											eq.append("*");
											eq.append(model.getCalculatedDescriptors().getColumn(k));
										}

										lvEquations.put("Latent variable " + (i + 1), eq.toString());
									}
									eData.supplementaryData.put("Equation", lvEquations);
								}
							}
						}

					}
					setFileName(model.name);
				}
			};

			wp.start();
			return redirect("longoperations/operationWaitingScreen.do?operation-id=" + wp.operationID);
		}
		else
		{
			Globals.setMarshallingOption(MarshallingOption.PROPERTY_UNITCATEGORY);
			Globals.setMarshallingOption(MarshallingOption.UNITCATEGORY_UNITS);
			return new WebModel(eData).setTemplate("export").getModelAndView();
		}
	}



	public ModelAndView attachment(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		// Get XML attachment of the model, for debugging purposes mostly /
		// Midnighter

		Model model = (Model) Globals.session().get(Model.class, getLongParam("id"));

		response.setContentType("application/xml");
		response.setHeader("Content-Disposition", "attachment; filename=\"model_" + model.name.replaceAll("\\s+", "_") + ".xml\"");

		OutputStream out = response.getOutputStream();

		if (assertParam("configuration"))
			out.write(model.attachment.getData());
		else
		{
			writeReadyModelXML(model, out);
		}
		out.flush();
		out.close();

		return null;
	}

	/**
	 * Export a ZIP with two files required for the standalone predictior: model configuration and model data XML files.
	 */
	public void exportStandaloneModel(HttpServletRequest request, HttpServletResponse response) throws IOException, JAXBException
	{

		if (Globals.userSession() == null || Globals.userSession().user == null || !Globals.userSession().user.isSuperUser())
			throw new UserFriendlyException("This feature is only available to the OCHEM administrator");

		Model model = getRequestedModel();

		if(model == null)throw new UserFriendlyException("The requested public model is not available in the database.");

		if(!model.isStandaloneExportable())throw new UserFriendlyException("The export of such types of models is currently not supported.");

		response.setContentType("application/zip");
		response.setHeader("Content-Disposition", "attachment; filename=\"standalone-model-" + model.publicId + ".zip\"");

		ZipOutputStream out = new ZipOutputStream(response.getOutputStream());

		out.putNextEntry(new ZipEntry("model-conf.ochem.xml"));
		Globals.createMarshaller(true).marshal(ExportableModel.create(model), out);
		out.closeEntry();

		out.putNextEntry(new ZipEntry("model-data.ochem.xml"));
		writeReadyModelXML(model, out);
		out.closeEntry();

		out.flush();
		out.close();
	}

	public void exportModelFile(HttpServletRequest request, HttpServletResponse response) throws IOException, JAXBException {
		if (Globals.userSession() == null || Globals.userSession().user == null || !Globals.userSession().user.isSuperUser())
			throw new UserFriendlyException("This feature is only available to the OCHEM administrator");

		Model model = getRequestedModel();

		if(model == null)throw new UserFriendlyException("The requested public model is not available in the database.");

		ModelAbstractConfiguration configuration =  ((CDSModelData)model.readyModelAttachment.getObject().modelData).methodSpecificData;

		if(configuration.saveModels() && configuration.getSavedModelAsBytes() != null) {

			response.setContentType("APPLICATION/OCTET-STREAM");
			response.setHeader("Content-Disposition", "attachment; filename=\"model-file." + model.publicId);
			OutputStream out = response.getOutputStream();
			out.write(configuration.getSavedModelAsBytes());
			out.flush();
			out.close();

		}else
			throw  new UserFriendlyException("no saved model is available");
	}


	private void writeReadyModelXML(Model model, OutputStream out) throws JAXBException
	{
		// Add AD assessment info to the exported model
		ReadyModelAttachment readyModelAttachment = model.readyModelAttachment.getObject();
		readyModelAttachment.adConfigurations = new ArrayList<ADConfiguration>();
		for (ModelMapping mm : model.modelMappings)
		{
			SetStatistics trainingSetStats = ModelStatistics.get(mm).sets.get(0);
			List<String> dmList = trainingSetStats.distancesToModel;
			if (!dmList.isEmpty())
			{
				ADConfiguration adConf = trainingSetStats.getADConfiguration(dmList.get(0), null, mm, null);
				adConf.epIds = null;
				readyModelAttachment.adConfigurations.add(adConf);
			}
		}

		Globals.createMarshaller(true).marshal(readyModelAttachment, out);
	}

	public ModelAndView exportModelXml(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		if (getParam("id").contains(","))
		{
			// Multiple models in a ZIP file
			String[] ids = getParam("id").split(",");
			response.setContentType("application/zip");
			response.setHeader("Content-Disposition", "attachment; filename=\""+ids.length+"_models_combined.ochem.zip\"");

			ZipOutputStream out = new ZipOutputStream(response.getOutputStream());

			for (int i = 0; i < ids.length; i++) 
			{
				Model model = (Model) Globals.session().get(Model.class, Long.valueOf(ids[i]));
				out.putNextEntry(new ZipEntry("" + (i+1) + "_model_" + model.name.replaceAll("\\s+", "_") + ".ochem.xml"));
				Marshaller jaxbMarshaller = Globals.jaxbContext.createMarshaller();
				jaxbMarshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
				jaxbMarshaller.marshal(ExportableModel.create(model), out);
				out.closeEntry();
			}

			out.flush();
			out.close();
		}
		else
		{
			// Single model
			Model model = (Model) Globals.session().get(Model.class, Long.valueOf(getLongParam("id")));

			response.setContentType("application/xml");
			response.setHeader("Content-Disposition", "attachment; filename=\"model_" + model.name.replaceAll("\\s+", "_") + ".ochem.xml\"");

			OutputStream out = response.getOutputStream();
			Marshaller jaxbMarshaller = Globals.jaxbContext.createMarshaller();
			jaxbMarshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
			jaxbMarshaller.marshal(ExportableModel.create(model), out);
			out.flush();
			out.close();
		}

		return null;
	}

	public ModelAndView generateGUID(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Model model = getRequestedModel();
		Privileges privileges = model.getPrivileges(request);

		if (!privileges.canView)
			throw new UserFriendlyException("You do not have privileges to generate GUID for this model");

		ModelIdentity identity = new ModelIdentity(model);
		Globals.session().save(identity);

		return new WebModel(identity).getModelAndView();
	}

	public ModelAndView approve(HttpServletRequest request, HttpServletResponse response) throws Exception {
		Model model = getRequestedModel();
		ModelOperation.approveModel(model, assertParam("publishedAndCited"), getIntParam("qualityGrade"));

		return new WebModel(new Alert("Thank you!\nThe model has been approved.")).getModelAndView();
	}

	public static PointSelector getPointSelector(HttpServletRequest request)
	{
		RequestParser parser = new RequestParser(request);
		PointSelector pointSelector = null;
		if (parser.assertParam("best-predictions"))
			pointSelector = new BestPredictionsPointSelector(1.0 * parser.getLongParam("best-predictions") / 100, true, 0);
		else if (parser.assertParam("dm-threshold"))
			pointSelector = new BestPredictionsPointSelector(Double.valueOf(parser.getParam("dm-threshold")));

		return pointSelector;
	}

	public ModelAndView getEstimatedStatistics(HttpServletRequest request, HttpServletResponse response) throws Exception {

		ModelProfileData profile = new ModelProfileData();
		ModelMapping mm = getRequestedModelMapping();
		ModelStatistics ms = ModelStatistics.get(mm);
		ms.validationSetId = "all";
		ms.actualizeStatistics(mm);
		ms.recalculateStatistics(mm);

		ApplicabilityDomain ad = new ApplicabilityDomain();
		ad.setModel(mm, null, null);

		for (int i = 0; i < mm.model.getValidationSets().size(); i++) {
			SetStatistics realStatistics = ms.getSetStatisticsByBasket(mm.model.getValidationSets().get(i).id);
			PredictionResults pr = new PredictionResults();
			pr.setPredictions(realStatistics, ad, mm, 0);
			pr.bootstrapStatistics();
			profile.addStats(mm.model.getValidationSets().get(i), realStatistics, pr.simulatedStatistics);
		}

		return new WebModel(profile).getModelAndView();
	}

	private ModelMapping getRequestedModelMapping() {
		Model model = getRequestedModel();
		if (model == null || assertParam("mapping_id"))
		{
			if (!assertParam("mapping_id"))
				throw new UserFriendlyException("The model you a trying to access is not available or has been deleted");
			return (ModelMapping) Globals.session().get(ModelMapping.class, getLongParam("mapping_id"));
		}

		return model.modelMappings.get(0);
	}

	private ModelApplier getModelApplier()
	{
		return (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
	}

	public static List<Long> getValidationSetIDs(HttpServletRequest request)
	{
		List<Long> ids = new ArrayList<Long>();
		for (int i = 0; i < 20; i++)
			if (request.getParameter("validationsetid" + i) != null && !"".equals(request.getParameter("validationsetid" + i)))
				ids.add(Long.valueOf(request.getParameter("validationsetid" + i)));
		if (!ids.isEmpty())
			return ids;
		else
			return null;

	}
}
