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

import java.io.Serializable;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.SessionVariable;
import qspr.dao.Repository;
import qspr.entities.Alert;import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.PendingTask;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.PropertyValue;
import qspr.entities.Unit;
import qspr.export.ExportableColumn;
import qspr.export.ExportableSet;
import qspr.export.ExportableSetConfiguration;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.frontend.modelapplier.ModelPrediction;
import qspr.frontend.modelapplier.ModelPredictionPiece;
import qspr.frontend.modelapplier.PredictionRow;
import qspr.interfaces.ProvidedConditions;
import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.metaserver.util.ShortCondition;
import qspr.workflow.utils.QSPRConstants;
import qspr.modelling.ApplicabilityDomain;
import qspr.modelling.CompoundsProvider;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.applier.ModelApplierTaskProcessor;
import qspr.modelling.applier.PredictionResults;
import qspr.modelling.configurations.ExternalCondition;
import qspr.util.ExportThread;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.useractions.EventFactory;
import com.eadmet.useractions.PredictionAction;
import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.NumericalValueStandardizer;

@Controller
@SuppressWarnings("unchecked")
public class ModelApplierController extends BrowserWrapper
{

	public ModelAndView apply(HttpServletRequest request, HttpServletResponse response)
	{
		if (assertParam("model"))
		{
			ModelApplier applier = new ModelApplier();
			Globals.setSessionAttribute(SessionVariable.MODEL_APPLIER, applier);
			String[] models = getParam("model").split(",");
			for (String modelId : models)
			{
				Model m = (Model) Globals.session().get(Model.class, Long.valueOf(modelId));
				applier.addModel(m);
			}
		}

		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		applier.reset();

		if (assertParam("scenario"))
			applier.scenario = PredictionScenario.PREDICTION_ONLY;

		if (assertParam("disable-cache"))
			applier.useCache = false;

		if (assertParam("force_cache")){
			applier.useCache = false;
			applier.forceUpdateDescriptorCache = true;
		}


		CompoundsProvider compoundsProvider = new CompoundsProvider();
		try
		{
			applier.compoundsProvider = compoundsProvider.parseUI();
		} catch (Exception e)
		{
			return new WebModel(new Alert(compoundsProvider.error)).setTemplate("model/providedata").getModelAndView();
		}

		if (applier.compoundsProvider.hasCompounds())
		{
			Repository.user.checkEligibility(compoundsProvider.getCompoundsNum(), QSPRConstants.APPLIER_BONUS);

			if (assertParam("preferred-server")) {
				applier.modelTasks.get(0).preferredServer = getParam("preferred-server").trim();
				applier.modelTasks.get(0).defaultTaskPriority = TaskPriority.EXTRA_HIGH;
			}

			if (applier.defaultConditions == null)
			{
				List<ExternalCondition> eDescs = applier.getNecessaryConditions();
				if (!eDescs.isEmpty())
				{
					applier.useCache = false; // no cache will be used if models are with external conditions; only predictions will be used
					Globals.setMarshallingOption(MarshallingOption.PROPERTY_UNITCATEGORY);
					Globals.setMarshallingOption(MarshallingOption.UNITCATEGORY_UNITS);
					Globals.setMarshallingOption(MarshallingOption.PROPERTY_OPTIONS_APPLIER);
					WebModel wm = new WebModel().addObjects(eDescs).setTemplate("model/provide-conditions");
					for (ShortCondition eDesc : eDescs)
						wm.listedObjects.add(Globals.session().get(Property.class, eDesc.id));
					return wm.getModelAndView();
				}
			}
	
			return new WebModel().setTemplate("model/apply").getModelAndView();
		}
		else
		{
			return new WebModel().setTemplate("model/providedata").getModelAndView();
		}
	}

	private void formModelPredicionPiece(ModelPredictionPiece pred, DataTable modelResult, ModelMapping mm, Integer index, List<ModelMapping> modelList, ModelApplierTaskProcessor mt) throws Exception
	{
		if ("error".equals(modelResult.getRow(index).status))
		{
			pred.value = "error";
			pred.error = modelResult.getRow(index).detailedStatus;
		}
		else
		{
			// Valid result - no error
			int colIndex = 0;
			int arrIndex = modelResult.getColumnIndex(QSPRConstants.INDIVIDUAL_PREDICTIONS);
			if (modelList.size() > 1)
			{
				colIndex = modelResult.getColumnIndex(QSPRConstants.PREDICTION_RESULT_COLUMN + mm._class);
				arrIndex = modelResult.getColumnIndex(QSPRConstants.INDIVIDUAL_PREDICTIONS + mm._class);
			}

			Double val = 0D;
			if (modelResult.getRowsSize() > 0)
				val = (Double) modelResult.getValue(index, colIndex);

			String cell = "";
			if (mm.property.isQualitative())
				cell += mm.model.attachment.getObject().getOptionFromPrediction(val, mm.property).name.toString();
			else
				cell += NumericalValueStandardizer.getSignificantDigitsStr(val, 2);

			if (mm.unit != null && mm.property.isNumeric())
				cell += " " + mm.unit.getName();

			ApplicabilityDomain ad = mt.getApplicabilityDomain(null, mm.getIndex());
			if (ad != null)
			{
				if (mm.property.isNumeric())
					cell += String.format(" &plusmn; %5.2f (%s = %5.2f, estimated RMSE = %5.2f)", ad.getPredictedError(index) * 1.96,
							ad.dmName, ad.getPredictedDM(index), ad.getPredictedError(index));
				else
					cell += " (" + Math.rint(ad.getPredictedError(index) * 100) + "% accuracy)";
				if (ad.dmThreshold != null && ad.getPredictedDM(index) > ad.dmThreshold)
					cell += "out_of_ad";
			}

			if (modelResult.getRow(index).getAttachment(QSPRConstants.CACHED) != null)
				cell += QSPRConstants.CACHED;

			pred.value = cell;

			if (arrIndex != -1 && modelResult.getRow(index).getValue(arrIndex) instanceof float[] && (float[])modelResult.getRow(index).getValue(arrIndex) != null) // TODO current fix for consensus network; second value should not be null
			{
				float[] ensemblePredictions = (float[])modelResult.getRow(index).getValue(arrIndex);
				pred.predictionVector = new ArrayList<Float>();
				for (float f : ensemblePredictions)
					pred.predictionVector.add(f);
			}
		}
	}

	public ModelAndView reslist(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		Map<Long, Integer[]> errorIndices = (Map<Long, Integer[]>)Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER_INDEXMAP);


		long listSize = 0;
		int pagenum = getPageNum();
		int pagesize = getPageSize(15);
		Integer[] mapList = null;

		if (assertParam("row_num"))
		{
			pagenum = 1;
			pagesize = 1;
		}

		if ("all".equals(request.getParameter("display-mode")) || request.getParameter("display-mode") == null)
		{
			listSize = applier.compoundsProvider.basket.getRowsSize();
		} else
			if ("errors".equals(request.getParameter("display-mode")))
			{
				mapList = errorIndices.get(-1L);
				listSize = mapList.length;
			}
			else
				if (request.getParameter("display-mode").startsWith("model"))
				{
					Long modelId = Long.valueOf(request.getParameter("display-mode").replaceAll("model-", ""));
					mapList = errorIndices.get(modelId);
					listSize = mapList.length;
				}

		@SuppressWarnings("rawtypes")
		List result = new ArrayList<PredictionRow>();

		applier.setResultOrdering(request.getParameter("ordering"), assertParam("asc"));

		for (int i = (pagenum - 1) * pagesize; i < Math.min(listSize, pagenum * pagesize); i++)
		{

			int index = 0;
			if (mapList == null)
				index = applier.order != null ? applier.order.get(i) : i;
			else
				index = mapList[i];

			if (assertParam("row_num"))
				index = getIntParam("row_num");

			PredictionRow row = new PredictionRow();

			ExperimentalProperty ep = applier.compoundsProvider.basket.entries.get(index).ep;// Observed
			if (ep.id != null && ep.property == null) // refetch
				ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, ep.id);

			if (ep.molecule != null && ep.molecule.id != null)
				row.id =  ep.molecule.id;
			else
				row.id = -1L;

			for (int j = 0; j < applier.modelTasks.size(); j++)
			{
				ModelApplierTaskProcessor mt = applier.modelTasks.get(j);
				mt.initialiseModel();
				Model m = mt.model;
				List<ModelMapping> modelList = m.modelMappings;

				ModelPrediction modelprediction = new ModelPrediction();
				modelprediction.modelId = mt.model.id;
				modelprediction.taskNum = j;
				modelprediction.trainingSetId = mt.model.trainingSet.id;

				if (mt.wndResult != null)
				{
					DataTable modelResult = mt.wndResult.ports.get(0);

					for (ModelMapping mm : modelList)
					{
						ModelPredictionPiece pred = new ModelPredictionPiece(mm.property.getName() + " (" + m.name + ") = ");
						pred.modelMappingId = mm.id;
						pred.tableRowNum = index;

						modelprediction.modelPredictionPieces.add(pred);

						formModelPredicionPiece(pred, modelResult, mm, index, modelList, mt);

						if (ep.id != null && ep.property != null && ep.property.id.equals(mm.property.id))
						{
							ModelPredictionPiece real = new ModelPredictionPiece(mm.property.getName() + "(measured) = ");
							real.value = ep.getStringValue();
							if (mm.property.isNumeric())
							{
								real.value += " " + ep.unit.getName();
								if (!ep.unit.equals(mm.unit))
									real.value += " = " + NumericalValueStandardizer.getSignificantDigits(ep.getConvertedValue(mm.unit)) + " " + mm.unit.getName();
							}
							modelprediction.modelPredictionPieces.add(real);
						}

					}
				}
				else
					modelprediction.error = mt.getStatus();

				row.modelPredictions.add(modelprediction);
			}

			result.add(row);
		}

		WebList list = new WebList();
		list.list = result;
		list.size = Long.valueOf(listSize).intValue();
		list.pageNum = pagenum;
		list.pageSize = pagesize;

		return new WebModel(list).getModelAndView();
	}

	public ModelAndView provideConditions(HttpServletRequest request, HttpServletResponse response) throws ParseException
	{
		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		if (assertParam("condition-id"))
		{
			ConditionSet defaultConditions = new ConditionSet();
			String[] condIds = request.getParameterValues("condition-id");
			String[] condValues = request.getParameterValues("condition-value");
			String[] condUnits = request.getParameterValues("condition-unit");

			for (int i = 0; i < condIds.length; i++)
			{
				PropertyValue pv = new PropertyValue();
				pv.property = Property.getById(Long.valueOf(condIds[i]));
				if (pv.property.isNumeric())
				{
					try{
						pv.value = new Double(condValues[i]);
					}catch(Exception e){ // trying to pass comma instead of point
						NumberFormat nf = NumberFormat.getInstance(Locale.GERMAN);
						pv.value = nf.parse(condValues[i]).doubleValue();
					}

					pv.unit = (Unit) Globals.session().get(Unit.class, Long.valueOf(condUnits[i]));
				}
				else
					pv.option = (PropertyOption) Globals.session().get(PropertyOption.class, Long.valueOf(condValues[i]));
				defaultConditions.values.add(pv);
			}
			applier.defaultConditions = defaultConditions;
			return new WebModel().setTemplate("model/apply").getModelAndView();
		}
		else
			return new WebModel().setTemplate("model/providedata").getModelAndView();
	}

	public ModelAndView getAppliedMoleculeId(HttpServletRequest request, HttpServletResponse response)
	{
		ModelApplier mbData = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		return new WebModel(mbData.compoundsProvider.basket.entries.get(getIntParam("num")).ep.molecule).getModelAndView();
	}

	public ModelAndView status(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		String status = "";
		if (applier == null)
			status = "Modeller is not initialized";
		else
		{
			StringBuilder statuses = new StringBuilder();
			for (ModelApplierTaskProcessor modelTask : applier.modelTasks)
			{
				modelTask.update();
				statuses.append(modelTask.getStatus() + "_$$_");
			}
			status = statuses.toString();

			if (applier.isReady())
				if (!applier.isError())
					status = "Finished";
				else
					status = applier.getErrorMessage();

		}
		return new WebModel(new Alert(status)).getModelAndView();
	}

	public ModelAndView start(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		for (ModelApplierTaskProcessor mTask : applier.modelTasks)
			mTask.initialiseModel();

		EventFactory.document("Prediction", new PredictionAction(applier), null);
		Globals.restartAllTransactions(true);
		applier.start();
		return new WebModel().getModelAndView();
	}

	public ModelAndView fetchlater(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Globals.setSessionAttribute(SessionVariable.MODEL_APPLIER, null); // Model applier will have to face the garbage collector
		System.gc();
		return new WebModel().setTemplate("model/fetch-later").getModelAndView();
	}

	public ModelAndView results(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ModelApplierSummary summary = new ModelApplierSummary();
		WebModel wm = new WebModel(summary);
		//		List<Object> objects = new ArrayList<Object>();

		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		Map<Long, Integer[]> indexMap = (Map<Long, Integer[]>)Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER_INDEXMAP);

		if (indexMap == null)
		{
			indexMap = new HashMap<Long, Integer[]>();
			Globals.setSessionAttribute(SessionVariable.MODEL_APPLIER_INDEXMAP, indexMap);
		}

		PendingTask task = null;

		if (assertParam("task"))
		{
			applier = new ModelApplier(task = (PendingTask) Globals.session().get(PendingTask.class, getLongParam("task")));
			Globals.setSessionAttribute(SessionVariable.MODEL_APPLIER, applier);
			indexMap = new HashMap<Long, Integer[]>();
			Globals.setSessionAttribute(SessionVariable.MODEL_APPLIER_INDEXMAP, indexMap);
		}

		List<Property> predictedBasketProperties = null;

		if (applier == null)
			return redirect("model/select.do");

		if (applier.compoundsProvider.basket.id != null)
			predictedBasketProperties = applier.compoundsProvider.basket.getProperty();

		List<Integer> allErrors = new ArrayList<Integer>();
		for (ModelApplierTaskProcessor modelTask : applier.modelTasks)
		{
			if (modelTask.wndResult != null && modelTask.wndResult.ports.get(0).getRowsSize() < applier.compoundsProvider.basket.getRowsSize())
				throw new UserFriendlyException("There was a problem applying your model. The number of returned compounds " + modelTask.wndResult.ports.get(0).getRowsSize() + " is not equal to the number of requested compounds " + applier.compoundsProvider.basket.getRowsSize());

			modelTask.initialiseModel();
			summary.models.add(modelTask.model);

			// Form lists of error indices
			List<Integer> localErrors = new ArrayList<Integer>();
			DataTable dt = modelTask.wndResult.ports.get(0);
			for (AbstractDataRow row : dt)
				if (row.isError())
				{
					localErrors.add(dt.currentRow);

					if (!allErrors.contains(dt.currentRow))
						allErrors.add(dt.currentRow);

					summary.thereAreErrors = true;
					modelTask.model.predictionErrors++;
				}
			indexMap.put(modelTask.model.id, localErrors.toArray(new Integer[0]));


			int propNum = 0;
			for (ModelMapping mm : modelTask.model.modelMappings)
			{
				ApplicabilityDomain applicabilityDomain = modelTask.getApplicabilityDomain(null, mm.getIndex());
				if (applicabilityDomain != null)
				{
					applicabilityDomain.modelMapping = (ModelMapping) Globals.session().get(ModelMapping.class, applicabilityDomain.modelMapping.id);
					if (applicabilityDomain.ssTraining.points.size() < 2000) // a temporary fix not to show AD plot for big models, see issue DEV-452
					{
						summary.ads.add(applicabilityDomain);
					}

					// If we have some predictions, estimate the average expected performance for the whole set
					if (modelTask.wndResult.ports.get(0).getRowsSize() >= 3)
					{
						PredictionResults pr = new PredictionResults();
						pr.setPredictions(modelTask.wndResult.ports.get(0), applicabilityDomain, mm, propNum);
						pr.bootstrapStatistics();
						summary.results.add(pr);
					}
				}

				if (predictedBasketProperties != null && predictedBasketProperties.contains(mm.property))
					wm.addParam("can-be-validation-set", modelTask.model.id.toString());

				propNum++;
			}

		}
		indexMap.put(-1L, allErrors.toArray(new Integer[0]));

		summary.conditions = applier.getConditionsDescription();
		if(summary.conditions == null && task != null)summary.conditions = task.description;
		if(summary.conditions == null) {
			Serializable conf = (Serializable) applier.modelTasks.get(0).model.attachment.getObject().configuration;
			if(conf instanceof ProvidedConditions && ((ProvidedConditions) conf).hasConditions()) {
				applier.defaultConditions = Repository.property.createDefaultConditions((ProvidedConditions)conf,false);
				summary.conditions = applier.getConditionsDescription();
				applier.defaultConditions = null;
			}
		}

		return wm.setTemplate("model/prediction-results").getModelAndView();
	}

	public ModelAndView exportPredictions(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		logger.info(MemoryUtils.memorySummary());

		ExportableSet eData = new ExportableSet();
		eData.uncheckColumn(ExportableColumn.DESCRIPTORS);
		eData.uncheckColumn(ExportableColumn.INTRODUCER);
		eData.uncheckColumn(ExportableColumn.MODIFIER);

		// To be able to convert "real" values here we need to obtain the list of involved properties to display them in UI
		// eData.properties = applier.getRealProperties() or something like that
		// Not highest priority right now

		final ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);

		if (assertParam("submit"))
		{
			logger.info(MemoryUtils.memorySummary());

			ExportThread eThread = applier.getExportThread(ExportableSetConfiguration.configureFromDialog(request), getParam("format"));
			eThread.start();

			return redirect("longoperations/operationWaitingScreen.do?operation-id=" + eThread.operationID);
		}
		else
		{
			// Make properties available from UI, export dialog
			Set<Property> props = new HashSet<Property>();
			for (ModelApplierTaskProcessor mTask : applier.modelTasks)
			{
				for (ModelMapping mm : mTask.model.modelMappings)
				{
					if (!Globals.session().contains(mm.property))
						mm.property = (Property) Globals.session().get(Property.class, mm.property.id);
					props.add(mm.property);
					mm.property.selectedUnit = mm.unit;
				}
			}
			Globals.setMarshallingOption(MarshallingOption.PROPERTY_UNITCATEGORY);
			Globals.setMarshallingOption(MarshallingOption.UNITCATEGORY_UNITS);
			eData.properties.addAll(props);
			return new WebModel(eData).setTemplate("export").getModelAndView();
		}
	}

	//Mixture-related code
	public ModelAndView provideMixture(HttpServletRequest request, HttpServletResponse response)
	{
		ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		applier.useCache = false;
		throw new UserFriendlyException("Not supported -- sorry");
		/*
		Property mixMol = Property.getByName(QSPRConstants.mixtureConditionMolId);
		Property mixFrac = Property.getByName(QSPRConstants.mixtureConditionMolarFraction);
		Molecule m = Repository.molecule.getMolecule(getLongParam("n-molecule"));

		Basket b = applier.compoundsProvider.getBasket();
		int i=0;
		for (BasketEntry be : b.entries)
		{
			i++;
			ConditionSet cs = null;

			if (be.ep.conditions == null)
				cs = new ConditionSet();
			else
				cs = (ConditionSet)Globals.session().get(ConditionSet.class, be.ep.conditions.id);

			if (cs.getValue(QSPRConstants.mixtureConditionMolId) == null)
			{
				cs.values.add(new PropertyValue(mixMol, "M" + m.mapping2.id));
				cs.values.add(new PropertyValue(mixFrac, Double.valueOf(getParam("n-molfrac"))));
			}
			Globals.session().evict(cs);
			if (i % 100 == 0)
				Globals.restartAllTransactions(true);
		}

		ExperimentalProperty ep = new ExperimentalProperty();
		ep.molecule = m;
		ConditionSet cs = new ConditionSet();
		cs.values.add(new PropertyValue(mixMol, "M0"));
		cs.values.add(new PropertyValue(mixFrac, 1D));
		ep.conditions = cs;
		applier.entries().add(new BasketEntry(ep));


		return new WebModel().setTemplate("model/apply").getModelAndView();
		 */

	}

	//Mixture-related code end


	@Override
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		// TODO Auto-generated method stub
		return null;
	}

	private static final Logger logger = LogManager.getLogger(ModelApplierController.class);

}

@XmlRootElement
class ModelApplierSummary
{
	@XmlAttribute
	String conditions;

	@XmlAttribute
	Boolean thereAreErrors;

	@XmlElement
	List<Model> models = new ArrayList<Model>();

	@XmlElement
	List<ApplicabilityDomain> ads = new ArrayList<ApplicabilityDomain>();

	@XmlElement
	List<PredictionResults> results = new ArrayList<PredictionResults>();	
}
