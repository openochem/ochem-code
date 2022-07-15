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
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.TimeoutException;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Hibernate;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.Transactional;
import qspr.dao.Repository;
import qspr.entities.Alert;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.PendingTask;
import qspr.entities.Property;
import qspr.frontend.WebModel;
import qspr.workflow.utils.QSPRConstants;
import qspr.modelling.ApplicabilityDomain;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.applier.ModelApplierTaskProcessor;
import qspr.util.AccessChecker;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.MAILERConstants;
import com.eadmet.utils.mailer.Mailer;

@Controller
public class PhyspropController extends ControllerWrapper
{
	private static transient final Logger logger = LogManager.getLogger(PhyspropController.class);

	/**
	 * Get the default set of models
	 */
	private List<Model> getDefaultModels()
	{
		List<Model> models = Repository.model.getFeaturedModels();

		// Sort alphabetically
		Collections.sort(models, new Comparator<Model>()
		{

			@Override
			public int compare(Model o1, Model o2)
			{
				return o1.featuredName.compareTo(o2.featuredName);
			}
		});

		return models;
	}

	/**
	 * Get a list of models given a comma-separated list of model IDs
	 */
	private List<Model> getModels(String ids)
	{
		List<Model> models = new ArrayList<Model>();
		String[] strIds = ids.split(",");
		for (String strId : strIds)
		{
			Model model = Repository.model.getById(Long.valueOf(strId));
			if (!models.contains(model))
				models.add(model);
		}
		return models;
	}

	private void forceLogin() {
		if (Globals.userSession() == null)
		{
			loginUser(null);
			Globals.userSession().license = true;
			Globals.userSession().disableQuota = true;
		}
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) 
	{
		forceLogin();
		return new WebModel().addObjects(getDefaultModels()).setTemplate("physprop").getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws IOException, TimeoutException {
		String action = getParam("action");

		if (action.equals("submit")){

			String data = getParam("moldata");
			logger.info(data);

			ModelApplier applier = new ModelApplier();
			applier.compoundsProvider.basket.addMolecule(data);

			for (Model model : getModels(getParam("modelIds"))) 
				applier.addModel(model);

			for (ModelApplierTaskProcessor modelTask : applier.modelTasks)
				modelTask.setUseCache();

			Globals.setSessionAttribute(SessionVariable.MODEL_APPLIER, applier);
		}

		return new WebModel().getModelAndView();
	}

	public ModelAndView start(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		forceLogin();
		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		if (applier == null){
			Mailer.notifyDevelopers("ePhyschem: model applier is null!", "Most probably the used browser does not store cookies.");
			throw new UserFriendlyException("It look like your browser does not support cookies or blocks our web-site. "
					+ "Please check your browser settings. If it does not help, please try reloading this page.\n If the error persists, please send a bug report to "+MAILERConstants.EMAIL_ADMIN);
		}
		applier.defaultTaskPriority = AccessChecker.getMaximumTaskPriority(Globals.userSession().user);
		applier.start();

		return new WebModel().getModelAndView();
	}

	public ModelAndView status(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		final ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		String status = "";

		if (applier == null)
			status = "Modeller is not initialized";
		else
		{
			if (applier.isReady())
			{
				status = "Finished";
				final DataTable[] dtRef = new DataTable[]{null};
				new Transactional()
				{

					@Override
					protected void wrapped() throws Exception
					{
						dtRef[0] = fillResultTable(applier);
					}
				}.execute();

				return new WebModel(dtRef[0]).getModelAndView();
			}

			int ready = 0;
			int applied = 0;
			for (ModelApplierTaskProcessor modelTask : applier.modelTasks)
			{
				//modelTask.model = (Model) Globals.session().get(Model.class, modelTask.model.id);
				modelTask.update();
				//status += modelTask.detailedStatus + "_$$_";

				applied = applied + modelTask.propertiesCount;
				if (modelTask.isReady())
					ready = ready + modelTask.propertiesCount;
			}

			status += ready + " / " + applied + " tasks are finished";

		}
		return new WebModel(new Alert(status)).getModelAndView();
	}

	List<String> results = new ArrayList<String>();

	private DataTable fillResultTable(ModelApplier applier) throws Exception
	{
		DataTable resTable = new DataTable();


		int index = 0;
		// create column headers
		for (int j = 0; j < applier.modelTasks.size(); j++)
		{
			ModelApplierTaskProcessor mt = applier.modelTasks.get(j);
			mt.initialiseModel();

			Hibernate.initialize(mt.model.modelMappings);
			for (ModelMapping mm : mt.model.modelMappings)
			{
				logger.info(mm.property.getName() + " (" + mt.model.name + ")");
				resTable.addColumn(mm.property.getName() + " (" + mt.model.name + ")");
			}
		}


		// value
		resTable.addRow();

		//fill columns with result values for each model
		int column = 0;
		for (int j = 0; j < applier.modelTasks.size(); j++)
		{
			ModelApplierTaskProcessor modelTask = applier.modelTasks.get(j);
			modelTask.initialiseModel();
			int lastColumn = column + modelTask.model.modelMappings.size();
			try {
				if (modelTask.isError())
					throw new UserFriendlyException(modelTask.getStatus());

				DataTable modelResult = modelTask.wndResult.ports.get(0);

				Model m = modelTask.model;
				List<ModelMapping> modelList = m.modelMappings;

				for (ModelMapping mm : modelList)
				{
					Hibernate.initialize(mm);
					String modelPrefix = mm.property.getName() + " (" + m.name + ") = ";

					if ("error".equals(modelResult.getRow(index).status))
					{
						resTable.setValue(column++, modelPrefix + modelResult.getRow(index).detailedStatus);
					}
					else
					{
						// Valid result - no error
						int colIndex = 0;
						if (modelList.size() > 1)
							colIndex = modelResult.getColumnIndex(QSPRConstants.PREDICTION_RESULT_COLUMN + mm._class);

						Double val = new Double(0);
						if (modelResult.getRowsSize() > 0)
							val = (Double) modelResult.getValue(index, colIndex);

						val = Math.rint(val * 100) / 100;
						String value = ""; 
						if (mm.property.isQualitative())
							value += mm.model.attachment.getObject().getOptionFromPrediction(val, mm.property).name.toString();
						else
							value += val.toString();

						String unit = "";
						if (mm.unit != null && mm.property.isNumeric())
							unit = mm.unit.getName();

						String ad = "";
						ApplicabilityDomain appd = modelTask.getApplicabilityDomain(null, mm.getIndex());
						if (appd != null)
							if (mm.property.isNumeric())
								ad += String.format(" &plusmn; %5.2f&#42;", appd.getPredictedError(index) * 0.97); // Midnighter on Dec 20, 2011 / 66% confidence interval, as discussed with Wolfram
							else
								ad += " (" + Math.rint(appd.getPredictedError(index) * 100) + "% accuracy)";

						resTable.setValue(column++, value + " " + unit + ad);

					}
				}
			} catch (Exception e) {
				e.printStackTrace();
				while (column < lastColumn)
				{
					resTable.setValue(column++, e.getMessage());
				}
			}
		}

		return resTable;
	}

	public ModelAndView results(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		WebModel wm = new WebModel();
		List<Object> objects = new ArrayList<Object>();
		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);

		if (assertParam("task"))
		{
			applier = new ModelApplier((PendingTask)Globals.session().get(PendingTask.class, getLongParam("task")));
			Globals.setSessionAttribute(SessionVariable.MODEL_APPLIER, applier);
		}

		List<Property> predictedBasketProperties = null;
		if (applier.compoundsProvider.basket.id != null)
			predictedBasketProperties = applier.compoundsProvider.basket.getProperty();

		for (ModelApplierTaskProcessor modelTask : applier.modelTasks)
		{
			modelTask.initialiseModel();
			objects.add(modelTask.model);
			for (ModelMapping mm : modelTask.model.modelMappings)
			{
				ApplicabilityDomain applicabilityDomain = modelTask.getApplicabilityDomain(null, mm.getIndex());
				if (applicabilityDomain != null)
				{
					applicabilityDomain.modelMapping = (ModelMapping) Globals.session().get(ModelMapping.class, applicabilityDomain.modelMapping.id);
					objects.add(applicabilityDomain);
				}

				if (predictedBasketProperties != null && predictedBasketProperties.contains(mm.property))
					wm.addParam("can-be-validation-set", modelTask.model.id.toString());
			}
		}

		return wm.setList(objects).getModelAndView();
	}
}
