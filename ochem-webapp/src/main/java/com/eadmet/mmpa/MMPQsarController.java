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

package com.eadmet.mmpa;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.SessionVariable;
import qspr.controllers.BrowserWrapper;
import qspr.controllers.NoFrame;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.PendingTask;
import qspr.entities.Property;
import qspr.entities.Session;
import qspr.export.CSVExportWriter;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.modelling.applier.ModelApplier;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.fingerprintindex.FingerprintIndex;
import com.eadmet.mmpa.domain.MMPAnnotationSet;
import com.eadmet.mmpa.domain.FuzzyValue;
import com.eadmet.mmpa.domain.MMPSubset;
import com.eadmet.mmpa.domain.MMPSubsetStatistics;
import com.eadmet.mmpa.domain.MMPTransformation;
import com.eadmet.mmpa.domain.MergedAnnotationSet;
import com.eadmet.mmpa.domain.MMPair;
import com.eadmet.mmpa.domain.ThinPair;
import com.eadmet.useractions.EventFactory;
import com.eadmet.utils.NumericalValueStandardizer;
import com.eadmet.utils.MAILERConstants;

@Controller
@SuppressWarnings("unchecked")
public class MMPQsarController extends BrowserWrapper
{

	public MMPQsarController() {
		sessionRequired = true;
	}

	public ModelAndView cleanCache(HttpServletRequest request, HttpServletResponse response) {
		MMPQueryService.clearQueryCache();
		Basket basket = Basket.getById(getLongParam("basket"));
		if(basket !=null)return basket(request, response);
		return model(request, response);
	}
	
	public ModelAndView model(HttpServletRequest request, HttpServletResponse response) {

		Property property = null;
		Model model = getRequestedModel();
		if (assertParam("property"))
			property = Property.getById(getLongParam("property"));

		Globals.setSessionAttribute(SessionVariable.MODEL_PROPERTY, property);

		MMPQuery mmp = new MMPQuery(model.trainingSet);
		WebModel wm = new WebModel(model);

		if(property != null){
			mmp.setProperty(property.id);
			wm.addObject(property);
		}

		return wm.addObject(MMPIndexingService.getInstance().getStatus(mmp)).
				setTemplate("matched-pairs/mmp-model").getModelAndView();
	}

	public ModelAndView basket(HttpServletRequest request, HttpServletResponse response) {
		Basket basket = Basket.getById(getLongParam("basket"));
		Property property = Property.getById(getLongParam("property"));
		return new WebModel(basket)
				.addObject(property)
				.addObject(MMPIndexingService.getInstance().getStatus(new MMPQuery(basket)))
				.setTemplate("matched-pairs/mmp-basket").getModelAndView();
	}

	public ModelAndView prediction(HttpServletRequest request, HttpServletResponse response) throws Exception {
		ModelMapping mm = getRequestedApplier().modelTasks.get(0).model.modelMappings.get(0);
		String dataset = getRequestedApplier().modelTasks.get(0).setDescription;
		return new WebModel(mm)
				.addParam("dataset", dataset)
				.setTemplate("matched-pairs/mmp-prediction").getModelAndView();
	}

	@NoFrame
	public ModelAndView getSignificanceData(HttpServletRequest request, HttpServletResponse response) throws IOException, ClassNotFoundException, Exception 
	{
		List<MMPTransformation>[] data 	 = 
				MMPSignificanceService.getSignificanceData(
						Basket.getBasket("Pyriformis Zhu complete"), 
						Property.getByName("log(IGC50-1)"), 
						new PendingTask[]{
								PendingTask.getById(376121L)
						});

		List<List<Object[]>> chart = new ArrayList<List<Object[]>>();
		for (List<MMPTransformation> list : data)
		{
			List<Object[]> rows = new ArrayList<Object[]>();

			for (MMPTransformation tr : list)
				rows.add(new Object[]{tr.id, Math.abs(tr.statistics.deltaMean), tr.statistics.pValue > 0 ? -Math.log(tr.statistics.pValue) / Math.log(10) : 16, tr.statistics});

			chart.add(rows);
		}

		return new ModelAndView("json", "object", chart);
	}

	@NoFrame
	public ModelAndView getSignificanceDataClass(HttpServletRequest request, HttpServletResponse response) throws IOException, ClassNotFoundException, Exception 
	{

		Basket b = Basket.getBasket("Aquatic Solubility (Classification)", Session.getLastSession(MAILERConstants.TESTER));
		Property p = Property.getByName("Aquatic Solubility Classification");
		PendingTask[] pt = new PendingTask[]{PendingTask.getById(413514)};

		List<MMPTransformation>[] data 	 = 	MMPSignificanceService.getSignificanceData(b, p, pt);

		List<List<Object[]>> chart = new ArrayList<List<Object[]>>();
		for (List<MMPTransformation> list : data)
		{
			if (list == null)
				continue;

			List<Object[]> rows = new ArrayList<Object[]>();

			for (MMPTransformation tr : list)
			{
				double np = tr.statistics.nNP * 1.0D / (tr.statistics.nNP + tr.statistics.nNN);
				double pn = tr.statistics.nPN * 1.0D / (tr.statistics.nPN + tr.statistics.nPP);
				double delta = 0;

				if (Double.isNaN(np))
					delta = pn;
				else if (Double.isNaN(pn))
					delta = np;
				else
					delta = Math.max(np, pn);

				rows.add(new Object[]{tr.id, delta, tr.statistics.pValue > 0 ? -Math.log(tr.statistics.pValue) / Math.log(10) : 16, tr.statistics});
			}

			chart.add(rows);
		}

		return new ModelAndView("json", "object", chart);
	}

	public ModelAndView getSignificanceChart(HttpServletRequest request, HttpServletResponse response) {
		return new WebModel()
				.setTemplate("matched-pairs/mmp-significance").getModelAndView();
	}

	@NoFrame
	public ModelAndView getChart(HttpServletRequest request, HttpServletResponse response) 
	{
		SubsetWithValues subset = getSubsetFromRequest();

		// Prepare the chart data
		List<ThinPair> pairs = subset.subset.getThinPairs(getRequestedTransformationId());
		logger.info("Preparing the combined DrDp chart data for " + pairs.size() + " pairs");

		List<Object[]> chartData = new ArrayList<Object[]>();

		Map<Integer, Integer> clusterMap = null;

		if (assertParam("clusters"))
		{
			clusterMap = new HashMap<Integer, Integer>();
			List<Object[]> clusters = FingerprintIndex.getInstance().getMMPClusters(pairs, 0.3);
			for (int i=0; i<clusters.size(); i++)
			{
				List<Integer> cl = (List<Integer>)clusters.get(i)[1];
				for (Integer mapping2 : cl)
					clusterMap.put(mapping2, i);
			}
		}

		for (int i = 0; i < pairs.size(); i++) 
		{
			ThinPair pair = pairs.get(i);
			if (!subset.valuesPredicted.containsKey(pair.mol1Id) || !subset.valuesPredicted.containsKey(pair.mol2Id))
				continue;
			String[] predictions = new String[]{
					NumericalValueStandardizer.getSignificantDigits(subset.valuesPredicted.get(pair.mol1Id).value), 
					NumericalValueStandardizer.getSignificantDigits(subset.valuesPredicted.get(pair.mol2Id).value)};
			double[] expValues = new double[]{subset.values.get(pair.mol1Id).value, subset.values.get(pair.mol2Id).value};
			chartData.add(new Object[]{pair.id, pair.mol1Id, pair.mol2Id, predictions, expValues, (clusterMap == null) ? 0 : clusterMap.get(pair.mol1Id)});
		}

		Map<String, Object> map = new HashMap<String, Object>();
		map.put("pairs", chartData);
		map.put("subset", subset.subset.id);

		return new ModelAndView("json", "object", map);
	}

	@NoFrame
	public ModelAndView getDeltaChartTransformations(HttpServletRequest request, HttpServletResponse response) throws Exception {
		MMPStatsRequest statsRequest = getTransformationsStatsRequest(request);

		Map<Long, MMPSubsetStatistics> statsMeasured = new HashMap<Long, MMPSubsetStatistics>();
		Map<Long, MMPSubsetStatistics> statsPredicted = new HashMap<Long, MMPSubsetStatistics>();

		List<MMPTransformation> trMeasured = MMPStatsService.getTransformationsStats(statsRequest);

		for (MMPTransformation mmpTransformation : trMeasured)
			statsMeasured.put(mmpTransformation.id, mmpTransformation.statistics);

		statsRequest.usePredictedValues = true;
		List<MMPTransformation> trPredicted = MMPStatsService.getTransformationsStats(statsRequest);
		for (MMPTransformation mmpTransformation : trPredicted)
			statsPredicted.put(mmpTransformation.id, mmpTransformation.statistics);

		List<Object[]> chartData = new ArrayList<Object[]>();
		for (Long trId : statsMeasured.keySet())
		{
			MMPSubsetStatistics ssMeasured = statsMeasured.get(trId);
			MMPSubsetStatistics ssPredicted = statsPredicted.get(trId);
			if (ssMeasured != null && ssPredicted != null)
				if (ssMeasured.pValue > 0)
					chartData.add(new Object[]{trId, ssMeasured.deltaMean, ssPredicted.deltaMean, -Math.log(ssMeasured.pValue) / Math.log(10)});
		}

		Map<String, Object> map = new HashMap<String, Object>();
		map.put("transformations", chartData);
		map.put("subset", statsRequest.subset.id);

		return new ModelAndView("json", "object", map);
	}


	public ModelAndView exportPairs(HttpServletRequest request, HttpServletResponse response) throws IOException {
		SubsetWithValues pairs = getSubsetFromRequest();

		OutputStream os = response.getOutputStream();

		response.setContentType("application/csv");
		response.setHeader("Content-Disposition", "attachment; filename=mmps.csv");
		CSVExportWriter eWriter = new CSVExportWriter();
		eWriter.os = os;
		eWriter.pw = new PrintWriter(os);
		eWriter.initialize();

		List<String> columns = new ArrayList<String>();
		columns.addAll(Arrays.asList(new String[]{"PAIR_ID", "TRANSFORMATION_ID", "MOL1", "MOL2", "SMIRKS"}));
		if (pairs.values != null)
		{
			columns.add("EXP_VALUE1");
			columns.add("EXP_VALUE2");
		}

		if (pairs.valuesPredicted != null)
		{
			columns.add("PREDICTION1");
			columns.add("PREDICTION2");
		}


		eWriter.writeRow(columns);

		for (ThinPair pair : pairs.subset.pairs)
		{
			List<Object> rows = new ArrayList<Object>();
			rows.add(pair.id);
			rows.add(pair.transformationId);
			rows.add("M" + pair.mol1Id);
			rows.add("M" + pair.mol2Id);
			rows.add(""+Globals.session().get(MMPTransformation.class,pair.transformationId)); // SMIRKS

			if (pairs.values != null)
			{
				rows.add(pairs.values.get(pair.mol1Id));
				rows.add(pairs.values.get(pair.mol2Id));
			}

			if (pairs.valuesPredicted != null)
			{
				rows.add(pairs.valuesPredicted.get(pair.mol1Id));
				rows.add(pairs.valuesPredicted.get(pair.mol2Id));
			}

			eWriter.writeRow(rows);
		}



		eWriter.flush();

		return null;
	}

	private SubsetWithValues getSubsetFromRequest() {
		SubsetWithValues result = new SubsetWithValues();

		MMPQuery query = new MMPQuery();
		MMPQueryService service = new MMPQueryService();

		if (getIntParam("similarity") != null && getIntParam("similarity") != 0)
			query.setSimilarity(getIntParam("similarity").shortValue());

		if (assertParam("basket")) {
			Basket basket = Basket.getById(getLongParam("basket"));
			Property property = Property.getById(getLongParam("property"));
			query.setBasket(basket.id).setProperty(property.id);
			result.subset = service.getPairsByBasket(query);
			result.values = MMPStatsService.getMolValues(basket.id, property.id, null);
		} else if (assertParam("prediction")) {
			query.applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
			//query.insideAD = true;//assertParam("insideAD");
			result.subset = service.getPairsFromPrediction(query);
			result.values = result.valuesPredicted = result.subset.molValues.values().iterator().next();
		}
		else if (getRequestedModel() != null)
		{
			Model model = getRequestedModel();
			Property property = getRequestedProperty();
			query.propertyId = property == null ? null: property.id;
			query.setBasket(model.trainingSet.id);
			result.subset = service.getPairsByBasket(query);
			result.valuesPredicted = MMPStatsService.getMolValues(model, true, property);
			result.values = MMPStatsService.getMolValues(model, false, property);
		}
		else
			throw new UserFriendlyException("Invalid request for value-delta chart");

		return result;
	}

	@NoFrame
	public ModelAndView getValueDeltaChart(HttpServletRequest request, HttpServletResponse response) throws Exception {
		SubsetWithValues res = getSubsetFromRequest();

		// Prepare the chart data
		logger.info("Preparing the combined DrDp chart data");
		List<ThinPair> pairs = res.subset.getThinPairs(getRequestedTransformationId());
		List<Object[]> chartData = new ArrayList<Object[]>();
		if (pairs != null)
			for (int i = 0; i < pairs.size(); i++) {
				ThinPair pair = pairs.get(i);
				if (!res.values.containsKey(pair.mol1Id) || !res.values.containsKey(pair.mol2Id))
					continue;
				double[] expValues = new double[]{res.values.get(pair.mol1Id).value, res.values.get(pair.mol2Id).value};
				chartData.add(new Object[]{pair.id, pair.mol1Id, pair.mol2Id, expValues});
				if (chartData.size() > 10000) // Do not overload UI
					break;
			}

		Map<String, Object> map = new HashMap<String, Object>();
		map.put("pairs", chartData);
		map.put("subset", res.subset.id);

		return new ModelAndView("json", "object", map);
	}

	/**
	 * Gather some statistics for the transformations with enough MMPs
	 */
	public ModelAndView getTransformationsStats(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		MMPStatsRequest statsRequest = getTransformationsStatsRequest(request);
		List<MMPTransformation> chosenTransformations = MMPStatsService.getTransformationsStats(statsRequest);

		EventFactory.document("MMP analysis", null, " has run MMP analysis ("+statsRequest+") and got " + chosenTransformations.size() + " significant transfomations");

		WebList wl = new WebList();
		wl.loadFromList(chosenTransformations, getPageNum(), getPageSize(5));

		return new BrowserModel().addFilter("subset", "" + statsRequest.subset.id, null).setObject(wl).getModelAndView();
	}

	public ModelAndView indexMolecules(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		MMPStatsRequest statsRequest = getTransformationsStatsRequest(request);
		MMPAnnotationService.indexMolecules(statsRequest.model != null? statsRequest.model.trainingSet : statsRequest.basket);
		return new WebModel().getModelAndView();
	}

	public ModelAndView annotate(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		MMPStatsRequest statsRequest = getTransformationsStatsRequest(request);
		String name = getParam("name");
		MMPAnnotationService.annotateRequestedTransformations(statsRequest, name);
		return new WebModel().getModelAndView();
	}

	//A special case for TransformationOptimizer merged set request
	public ModelAndView listAnnotationSetsMerged(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		String[] sets = request.getParameterValues("set");
		List<Long> setIds = new ArrayList<Long>();
		for (String set : sets)
			setIds.add(Long.valueOf(set));
		MergedAnnotationSet mas = new MergedAnnotationSet(setIds);
		MMPAnnotationService.calculateAnnotationCounts(mas);
		return new WebModel(mas).getModelAndView();
	}

	public ModelAndView listAnnotationSets(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		boolean published = "1".equals(getParam("published"));
		Criteria c = Globals.session().createCriteria(MMPAnnotationSet.class);
		Disjunction d = Restrictions.disjunction();
		d.add(Restrictions.eq("user", Globals.userSession().user));
		if (published)
			d.add(Restrictions.eq("featured", true));
		c.add(d);
		List<MMPAnnotationSet> l = c.list();

		if (assertParam("includeDefault"))
		{
			MMPAnnotationSet def = new MMPAnnotationSet();
			def.name = "Default";
			def.id = 0L;
			Globals.session().evict(def);
			l.add(0, def);
		}

		if (assertParam("counts"))
		{
			for (MMPAnnotationSet set : l) 
				MMPAnnotationService.calculateAnnotationCounts(set);
		}

		WebList wl = new WebList().loadFromList(l, 1, 100);
		return new WebModel(wl).getModelAndView();
	}

	public ModelAndView annotationSetProfile(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Long asId = getLongParam("id");
		MMPAnnotationSet as = (MMPAnnotationSet)Globals.session().get(MMPAnnotationSet.class, asId);

		if (as == null)
			throw new UserFriendlyException("MMP annotation set not found");

		boolean notMySet = Globals.userSession() == null || Globals.userSession().user == null || !as.user.equals(Globals.userSession().user);
		if (notMySet && !as.published)
			throw new UserFriendlyException("You cannot to this private set of transformations created by " + as.user.login);

		if (assertParam("action"))
		{
			if (notMySet)
				throw new UserFriendlyException("You are not authorized to modify this set of transformations created by " + as.user.login);

			String action = getParam("action");
			if ("rename".equals(action))
			{
				as.name = getParam("name");
				Globals.session().saveOrUpdate(as);
			} else
				if ("delete".equals(action))
				{
					Globals.session().createQuery("delete MMPTransformationAnnotation where annotationSet=:as").setParameter("as", as).executeUpdate();
					Globals.session().delete(as);
				} else
					if ("propdelete".equals(action))
					{
						Long propId = getLongParam("property");
						Property p = (Property)Globals.session().get(Property.class, propId);
						Globals.session().createQuery("delete MMPTransformationAnnotation where annotationSet=:as and property=:prop").setParameter("as", as).setParameter("prop", p).executeUpdate();
					}
		} else
		{
			MMPAnnotationService.calculateAnnotationCounts(as);
		}
		return new WebModel(as).setTemplate("matched-pairs/annotation-set-profile").getModelAndView();
	}

	@NoFrame
	public ModelAndView getFragmentGraph(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		MMPStatsRequest statsRequest = getTransformationsStatsRequest(request);
		List<MMPTransformation> chosenTransformations = MMPStatsService.getTransformationsStats(statsRequest);

		for (MMPTransformation transformation : chosenTransformations)
			if (transformation.statistics.deltaMean < 0)
				transformation.inversed = true;

		List<long[]> edges = MMPFragmentService.getFragmentGraph(chosenTransformations, getIntParam("cluster"), getLongParam("fragment"));
		return new ModelAndView("json", "object", edges);
	}

	public MMPStatsRequest getTransformationsStatsRequest(HttpServletRequest request) throws Exception {
		Globals.setMarshallingOption(MarshallingOption.PAIR_TRANSFORMATION);

		MMPStatsRequest statsRequest = new MMPStatsRequest();
		//statsRequest.maxPairs = 300000L;
		statsRequest.maxAtoms = getIntParam("maxAtoms");

		if (assertParam("prediction"))
		{
			statsRequest.applier = getRequestedApplier();
			statsRequest.bootstrapStatistics = assertParam("bootstrapStatistics");
			//statsRequest.insideAD = true;//assertParam("insideAD");
		}
		else if (assertParam("basket")) {
			statsRequest.basket = Basket.getById(getLongParam("basket"));
			statsRequest.property = Property.getById(getLongParam("property"));
		} else
		{
			statsRequest.model = getRequestedModel();	
		}

		if (assertParam("property"))
			statsRequest.property = getRequestedProperty();

		if (getIntParam("similarity") != null && getIntParam("similarity") > 0)
			statsRequest.similarity = getIntParam("similarity").shortValue();

		if (assertParam("minPairs"))
			statsRequest.minPairs = getIntParam("minPairs");

		if (assertParam("transformationPattern") && getParam("transformationPattern").length()>0)
			statsRequest.transformationPattern = getParam("transformationPattern");

		if (assertParam("pValue"))
			statsRequest.pValue = getDoubleParam("pValue");

		return statsRequest;
	}

	private Long getRequestedTransformationId() {
		if (assertParam("transformation"))
			return getLongParam("transformation");
		else if (assertParam("transformationPair"))
			return MMPair.get(getLongParam("transformationPair")).transformation.id;
		else
			return null;
	}

	public Property getRequestedProperty(){
		if (!assertParam("property"))
			return (Property) Globals.getSessionAttribute(SessionVariable.MODEL_PROPERTY);
		return Property.getById(getLongParam("property"));
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

	private ModelApplier getRequestedApplier() throws Exception {
		long taskId = getLongParam("prediction");
		if (taskId == 1)
			return (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
		else
			return ModelApplier.getFromTask(taskId);
	}

	private static final Logger logger = LogManager.getLogger(MMPQsarController.class);

	@Override
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		// TODO Auto-generated method stub
		return null;
	}
}

class SubsetWithValues {
	MMPSubset subset;
	Map<Integer, FuzzyValue> values;
	Map<Integer, FuzzyValue> valuesPredicted;
}
