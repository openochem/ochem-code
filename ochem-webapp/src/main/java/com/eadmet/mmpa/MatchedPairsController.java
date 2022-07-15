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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.SessionVariable;
import qspr.ThreadScope;
import qspr.controllers.BrowserWrapper;
import qspr.controllers.NoFrame;
import qspr.dao.Repository;
import qspr.entities.Mapping2;
import qspr.entities.Mapping2.MMPAIndex;
import qspr.entities.Property;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.modelling.applier.ModelApplier;
import qspr.util.MoleculePeer;
import qspr.util.RequestParser;
import qspr.util.ValueDistribution;
import qspr.util.WrapperThread;

import com.eadmet.fingerprintindex.FingerprintIndex;
import com.eadmet.mmpa.domain.MMPAnnotationSet;
import com.eadmet.mmpa.domain.MMPFragment;
import com.eadmet.mmpa.domain.MMPSubset;
import com.eadmet.mmpa.domain.MMPTransformation;
import com.eadmet.mmpa.domain.MMPair;
import com.eadmet.mmpa.domain.ThinPair;
import com.eadmet.mmpa.domain.MMPTransformationAnnotation;
import com.eadmet.utils.PageInfo;

@Controller
@SuppressWarnings("rawtypes")
public class MatchedPairsController extends BrowserWrapper
{
	private static transient final Logger logger = LogManager.getLogger(MatchedPairsController.class);

	public MatchedPairsController() {
		sessionRequired = true;
	}

	public ModelAndView transformations(HttpServletRequest req, HttpServletResponse res) {
		return new WebModel().setTemplate("matched-pairs/transformations").getModelAndView();
	}

	public ModelAndView clearCache(HttpServletRequest req, HttpServletResponse res) {
		MMPQueryService.clearQueryCache();
		MMPFragment.clearCache();
		MMPAnnotationService.countsCache.clear();
		MMPTransformation.clearCache();
		return redirect("matchedpairs/status.do");
	}

	/**
	 * Get the list of all or part of MMP transformations
	 * @throws Exception 
	 */
	public ModelAndView transformationsList(HttpServletRequest req, HttpServletResponse res) throws Exception {

		ThreadScope.setStatus("Loading the list of molecular transformations..");
		TransformationFilter filter = getFilter(req);

		WebList wl = new WebList().useEntity(MMPTransformation.class);
		MMPSubset subset = getQuerySubset();
		if (subset != null)
		{
			wl = getPaginatedTransformations(subset);
			return new BrowserModel().addFilter("subset", "" + subset.id, null).setObject(wl).getModelAndView();
		}
		else {
			MMPQueryService s = new MMPQueryService();
			filter.pageInfo = new PageInfo(getPageNum(), getPageSize(50));
			List<MMPTransformation> l = s.getTransformationsPage(filter);
			wl.loadFromList(l);
			wl.size = filter.pageInfo.size;
			wl.pageNum = filter.pageInfo.pageNum;
			wl.pageSize = filter.pageInfo.pageSize;
		}

		if (assertParam("getPairs"))
			Globals.removeMarshallingOption(MarshallingOption.PAIR_TRANSFORMATION);

		return new WebModel(wl).getModelAndView();
	}

	public ModelAndView getTransformation(HttpServletRequest req, HttpServletResponse res) {
		MMPTransformation transformation = MMPTransformation.get(getRequestedTransformationId());
		return new WebModel(transformation).getModelAndView();

	}

	public ModelAndView annotateProperty(HttpServletRequest req, HttpServletResponse res) {
		final long propertyId = getLongParam("id");
		WrapperThread thread = new WrapperThread()
		{
			@Override
			public void wrapped() throws Exception
			{
				MMPAnnotationService.annotateProperty(Property.getById(propertyId));
			}
		};

		thread.newOperation().start();

		return redirect("longoperations/operationWaitingScreen.do?operation-id=" + thread.operation.operationId);
	}

	public ModelAndView scheduleForIndexation(HttpServletRequest req, HttpServletResponse res) {
		MMPQuery query = createQuery();
		MMPIndexingService.getInstance().scheduleForIndexation(query);
		return new WebModel().getModelAndView();
	}

	/**
	 * Get a subset of MMPs relevant based on the query URL
	 * (for a basket, for a property or for a model applier)
	 * @throws  
	 */
	private MMPSubset getQuerySubset() 
	{
		MMPQuery query = new MMPQuery();

		if (getIntParam("similarity") != null && getIntParam("similarity") > 0)
			query.setSimilarity(getIntParam("similarity").shortValue());

		if (assertParam("subset"))
			return new MMPQueryService().getCachedSubset(getParam("subset"));
		else if (assertParam("basket")) 
		{
			query.setBasket(getLongParam("basket"));
			return new MMPQueryService().getPairsByBasket(query);
		}
		else if (assertParam("property"))
		{
			query.setProperty(getLongParam("property"));
			query.affectedPairsOnly = assertParam("affectedPairs");
			query.transformationId = getLongParam("transformation");
			query.requestData = query.transformationId != null;
			query.significantTransformationsOnly = assertParam("significantTransformations");

			if (assertParam("onlyIncreasing"))
				query.propertyChangeDirection = 1;
			else if (assertParam("onlyDecreasing"))
				query.propertyChangeDirection = -1;

			return new MMPQueryService().getPairsByProperty(query);
		}
		else if (assertParam("predictor"))
			return new MMPQueryService().getPairsFromPrediction(new MMPQuery((ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER), assertParam("insideAD")));
		else 
			return null;
	}

	private MMPQuery createQuery() {
		MMPQuery query = new MMPQuery();
		if (assertParam("basket"))
			query.basketId = getLongParam("basket");
		if (assertParam("tag"))
			query.tagId = getLongParam("tag");
		if (assertParam("property"))
			query.propertyId = getLongParam("property");
		if (!query.hasFilters())
			return null;
		return query;
	}

	@SuppressWarnings("unchecked")
	private WebList getPaginatedTransformations(MMPSubset subset) {
		WebList wl = new WebList();
		Long[] trIds = subset.getTransformationIds().toArray(new Long[]{});
		wl.loadFromList(Arrays.asList(trIds), getPageNum(), getPageSize(10));
		List transformations = new ArrayList<MMPTransformation>();
		for (int i = 0; i < wl.list.size(); i++)
		{
			MMPTransformation transformation = (MMPTransformation.get((Long) wl.list.get(i)));
			transformation.filteredPairs = subset.getPairs(transformation.id);
			transformations.add(transformation);
		}

		wl.list = transformations;

		return wl;
	}

	@SuppressWarnings("unchecked")
	private WebList getPaginatedPairs(MMPSubset subset) {
		WebList wl = new WebList();

		Long transformationId = null;
		if (assertParam("transformation"))
			transformationId = getLongParam("transformation");
		else if (assertParam("transformationPair"))
			transformationId = MMPair.get(getLongParam("transformationPair")).transformation.id;

		List l = subset.getThinPairs(transformationId);

		if (l == null)
			l = new ArrayList();

		wl.loadFromList(l, getPageNum(), getPageSize(100));

		List pairs = new ArrayList<MMPair>();

		for (int i = 0; i < wl.list.size(); i++)
			pairs.add(((ThinPair) wl.list.get(i)).getPair());

		wl.list = pairs;
		subset.fillPairsWithData(pairs);

		return wl;
	}

	//TODO change to automatic
	public ModelAndView transformationProfile(HttpServletRequest req, HttpServletResponse res) {
		String profile = "transformationProfile";

		String what = canExecuteLongRequest(profile, 1);

		if(what != null)
			logger.info("starting: " + what);
		else
		{
			logger.info("cannot start: " + profile);
			return redirect("user/pleaseregister.do");
		}

		Globals.setMarshallingOption(MarshallingOption.TRANSFOMATION_ANNOTATIONS);
		MMPTransformation transformation = MMPTransformation.get(getLongParam("id"));
		ModelAndView v = new WebModel(transformation).setTemplate("matched-pairs/transformation-profile").getModelAndView();
		cleanRequest(profile);
		logger.info("cleaning: " + profile);
		return v;
	}

	public ModelAndView getPairPredictions(HttpServletRequest req, HttpServletResponse res) {
		// TODO
		//MMPSubset subset = new MMPQueryService().getPairsFromPrediction(applier)
		return null;
	}

	public ModelAndView getClusters(HttpServletRequest req, HttpServletResponse res)
	{
		MMPSubset subset = getQuerySubset();
		List<ThinPair> pairs = subset.getThinPairs(getRequestedTransformationId());
		List<Object[]> clusters = FingerprintIndex.getInstance().getMMPClusters(pairs, 0.3);
		return new ModelAndView("json", "object", clusters);
	}

	@NoFrame
	public ModelAndView getHistogram(HttpServletRequest req, HttpServletResponse res) {
		List<Double> deltas;
		if (assertParam("subset"))
			deltas = MMPStatsService.getPairDeltas(getParam("subset"), getLongParam("transformation"));
		else
			deltas = MMPStatsService.getPairDeltas(getLongParam("property"), getLongParam("transformation"));
		ValueDistribution vd = new ValueDistribution();
		vd.setValues(deltas);
		return new ModelAndView("json", "object", vd);
	}

	public ModelAndView getPairs(HttpServletRequest req, HttpServletResponse res) {
		String profile = "getPairs";

		if(canExecuteLongRequest(profile, 1) == null)return redirect("user/pleaseregister.do");

		WebList wl = new WebList();

		MMPSubset subset = getQuerySubset();
		if (subset != null)
			wl = getPaginatedPairs(subset);
		else
		{
			setStatus("Loading pairs...");
			Criteria c = Globals.session().createCriteria(MMPair.class);
			if (assertParam("transformation"))
				c.add(Restrictions.eq("transformation", MMPTransformation.get(getLongParam("transformation"))));

			if (assertParam("molecule"))
			{
				Globals.setMarshallingOption(MarshallingOption.PAIR_TRANSFORMATION);
				Mapping2 mp2 = Repository.molecule.getMapping2(getIntParam("molecule"));
				c.add(Restrictions.or(Restrictions.eq("molecule1", mp2), Restrictions.eq("molecule2", mp2)));
			}

			if (assertParam("similarity"))
			{
				Integer s = getIntParam("similarity");
				if (s > 0)
					c.add(Restrictions.gt("similarity", s.shortValue()));
			}

			wl.loadFromCriteria(c, getPageNum(), getPageSize(100));
		}

		ModelAndView v = new WebModel(wl).getModelAndView();
		cleanRequest(profile);
		return v;
	}



	public ModelAndView status(HttpServletRequest req, HttpServletResponse res) {

		MMPIndexingService service = MMPIndexingService.getInstance();
		if (assertParam("recalculateCounts"))
		{
			service.recalculatePairCounts();
			Globals.restartAllTransactions(true);
		}

		if (assertParam("resubmit"))
		{
			Globals.session().createSQLQuery("update Mapping2 set mmpaIndexStatus=1 where mmpaIndexStatus=:status")
			.setParameter("status", "stuck".equals(getParam("resubmit")) ? MMPAIndex.MMP_SUBMITTED.ordinal() : MMPAIndex.MMP_FAILED.ordinal())
			.executeUpdate();
			Globals.restartAllTransactions(true);
			return redirect("matchedpairs/status.do");
		}

		if (assertParam("deleteInfrequent"))
		{
			service.recalculatePairCounts();
			Globals.restartAllTransactions(true);

			service.deleteInfrequentTransformations();
			Globals.restartAllTransactions(true);
		}

		if (assertParam("deleteDublicatePairs"))
		{
			service.deleteDublicatedPairs();
			service.recalculatePairCounts();
		}

		MMPStatus status = MMPIndexingService.getInstance().getStatus(createQuery());
		return new WebModel(status).setTemplate("matched-pairs/status").getModelAndView();
	}

	public ModelAndView startFragmentation(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		new WrapperThread()
		{

			@Override
			public void wrapped() throws Exception
			{
				//MMPIndexingService.getInstance().submitIndexingTasks();
			}
		}.start();
		return status(req, res);
	}

	// Copy-paste method. Refactor this piece
	private Long getRequestedTransformationId() {
		if (assertParam("transformation"))
			return getLongParam("transformation");
		else if (assertParam("transformationPair"))
			return MMPair.get(getLongParam("transformationPair")).transformation.id;
		else
			return null;
	}

	public ModelAndView startIndexing(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		new WrapperThread()
		{
			@Override
			public void wrapped() throws Exception
			{
				MMPIndexingService.getInstance().fetchIndexingTasks();
			}
		}.start();
		return status(req, res);
	}

	@NoFrame
	public ModelAndView getFragmentGraph(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		List<long[]> edges;
		/*if (assertParam("mmp-optimiser"))
		{
			MMPOptimiserProcessor optimiser = ((MMPOptimiserProcessor)Globals.getSessionAttribute(SessionVariable.MMPOPTIMISER));
			edges = MMPFragmentService.getFragmentGraph(optimiser.getAppliedTransformations(assertParam("effective-only")), getIntParam("cluster"), getLongParam("fragment"));
		}
		else */
		if (assertParam("explicitFilter"))
		{
			TransformationFilter filter = getFilter(request);
			edges = MMPFragmentService.getFragmentGraph(filter, getIntParam("cluster"), getLongParam("fragment"));
		}
		else if (assertParam("property") && !assertParam("basket"))
			edges = MMPFragmentService.getFragmentGraph(getLongParam("property"), getIntParam("cluster"), getLongParam("fragment"));
		else
			return new MMPQsarController().getFragmentGraph(request, response);
		return new ModelAndView("json", "object", edges);
	}


	public ModelAndView fragmentGraph(HttpServletRequest req, HttpServletResponse res) {
		return new WebModel().setTemplate("matched-pairs/fragment-graph").getModelAndView();
	}

	@Override
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		// TODO Auto-generated method stub
		return null;
	}

	public static  TransformationFilter getFilter(HttpServletRequest request) throws Exception {
		RequestParser parser = new RequestParser(request);
		TransformationFilter filter = new TransformationFilter();
		if (parser.assertParam("minCount"))
			filter.minCount = parser.getLongParam("minCount");
		if (parser.assertParam("exactCount"))
			filter.exactCount = parser.getLongParam("exactCount");

		//filter.fragId = 3L;
		filter.splitSets = "true".equals(parser.getParam("split-sets"));

		if (parser.assertParam("mol"))
		{
			Mapping2 mol = MoleculePeer.fetchFromString(parser.getParam("mol"), null).mapping2;
			if (mol.mmpaIndexStatus != MMPAIndex.MMP_INDEXED)
				MMPIndexingService.getInstance().instantIndex(mol);
			//throw new UserFriendlyException("The fragments for the molecule haven not yet been indexed (indexing status="+mol.mmpaIndexStatus+"). Please, try again later");
			filter.moleculeId = mol.id;
		}

		for (Object okey : request.getParameterMap().keySet())
		{
			String key = (String) okey;
			if (key.startsWith("p-"))
			{
				String[] parts = key.split("-");
				String propertyId = parts[1];
				String[] setIds = parts[3].split("\\.");
				String effect = parts[4];
				for (String setId : setIds) 
				{
					boolean[] effects = filter.getEffects(Long.valueOf(setId), Long.valueOf(propertyId));
					if ("increase".equals(effect))
						effects[2] = true;
					if ("decrease".equals(effect))
						effects[0] = true;
					if ("noeffect".equals(effect))
						effects[1] = true;
				}
			}
		}

		filter.clean();
		return filter;
	}

	public ModelAndView annotatedTransformations(HttpServletRequest req, HttpServletResponse res) {
		Property property = Property.getById(getLongParam("property"));
		MMPAnnotationSet set = MMPAnnotationSet.getById(getLongParam("set"));

		return new WebModel()
				.setTemplate("matched-pairs/annotated-transformations")
				.addObject(set)
				.addObject(property)
				.getModelAndView();
	}

	public ModelAndView annotatedTransformationsList(HttpServletRequest req, HttpServletResponse res) {

		Globals.setMarshallingOption(MarshallingOption.PAIR_TRANSFORMATION);

		Property property = Property.getById(getLongParam("property"));
		MMPAnnotationSet set = MMPAnnotationSet.getById(getLongParam("set"));

		Criteria c = Globals.session().createCriteria(MMPTransformationAnnotation.class)
				.add(Restrictions.eq("annotationSet", set))
				.add(Restrictions.eq("property", property));

		WebList wl = new WebList();
		wl.loadFromCriteria(c, getPageNum(), getPageSize(10));

		return new BrowserModel()
				.setObject(wl)
				.getModelAndView();
	}
}
