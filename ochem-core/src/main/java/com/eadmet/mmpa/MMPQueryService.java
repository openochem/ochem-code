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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.Query;
import org.hibernate.criterion.Criterion;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Junction;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.exception.GenericJDBCException;
import org.hibernate.sql.JoinType;
import org.hibernate.type.IntegerType;
import org.hibernate.type.LongType;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.ThreadScope;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.Mapping2Filter;
import qspr.entities.ModelMapping;
import qspr.entities.Property;
import qspr.modelling.ApplicabilityDomain;
import qspr.modelling.ModelStatistics;
import qspr.modelling.applier.ModelApplierTaskProcessor;
import qspr.modelling.applier.PredictionResults;
import qspr.modelling.applier.PropertyPrediction;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.mmpa.domain.MMPAnnotationSet;
import com.eadmet.mmpa.domain.FuzzyValue;
import com.eadmet.mmpa.domain.MMPSubset;
import com.eadmet.mmpa.domain.MMPTransformation;
import com.eadmet.utils.MemoryUtils;

/**
 * A service for querying MMPs according to different criteria (basket, property, ..)\
 * 
 * @author midnighter
 *
 */
@SuppressWarnings("unchecked")
public class MMPQueryService extends Thread
{
	/**
	 * A cache: <query key> -> [<transformation> -> <pairs>]
	 */
	private static Map<String, MMPSubset> queryCache = new HashMap<String, MMPSubset>();

	static Integer count = null;

	static public void clearQueryCache()
	{
		queryCache.clear();
	}

	public MMPSubset getPairs(MMPQuery query) 
	{
		if (query.basketId != null)
			return getPairsByBasket(query);
		else if (query.propertyId != null)
			return getPairsByProperty(query);
		else if (query.model != null)
			return getPairsFromPredictedBasket(query);
		else
			throw new UserFriendlyException("Invalid MMP query");
	}

	public MMPSubset getPairsByBasket(Long basketId) 
	{
		MMPQuery query = new MMPQuery();
		query.basketId = basketId;
		return getPairsByBasket(query);
	}

	public MMPSubset getPairsByBasket(MMPQuery query) {
		String key = "b-" + query.basketId;
		if (query.propertyId != null)
			key += "-p-" + query.propertyId;
		if (query.similarity != null)
			key += "-sim-" + query.similarity;

//		queryCache.clear();
		
		if (queryCache.containsKey(key))
			return queryCache.get(key);

		MMPSubset mmps = new MMPSubset();

		query.filterId = createMoleculeFilter(query);
		mmps = queryByFilter(query);
		Mapping2Filter.clear(query.filterId);

		addToCache(key, mmps);

		return mmps;
	}

	public int createMoleculeFilter(MMPQuery query) {
		int filterID = Mapping2Filter.generateFilterID();

		logger.info("Creating the molecule filter");
		String tables = "ExperimentalProperty join Molecule using (molecule_id)";
		List<String> where = new ArrayList<String>();
		if (query.propertyId != null)
			where.add("property_id = " + query.propertyId);
		if (query.basketId != null)
		{
			tables += " join BasketEntry using (exp_property_id)";
			where.add("basket_id=" + query.basketId);
		}

		String queryStr = "insert into Mapping2Filter(filter_id, mapping2_id) (select "+filterID+", mapping2_id from " + tables + " where " + StringUtils.join(where, " and ");

		if (query.basketId == null)
			queryStr += " and rights=2 and exp_property_id=first_entry";

		queryStr += "  group by mapping2_id)";

		Query q = Globals.session().createSQLQuery(queryStr);

		try
		{
			q.executeUpdate();
		} catch (GenericJDBCException e)
		{
			if (e.getSQLException().getMessage().contains("is full"))
			{
				logger.warn("The filter table is out of capacity. Clearing the table and retrying...");
				Globals.session().createSQLQuery("delete from Mapping2Filter").executeUpdate();
				Globals.restartAllTransactions(true);
				Globals.session().createSQLQuery(queryStr).executeUpdate();
			}
		}

		return filterID;
	}

	public MMPSubset queryByFilter(MMPQuery mmpQuery) {
		MMPSubset mmps = new MMPSubset();

		//int filterID = createMoleculeFilter(mmpQuery);

		String queryStr = "select transformation_id, mmp_id, mol1, mol2 from MMPair ";

		if (mmpQuery.significantTransformationsOnly)
			queryStr += " natural join MMPTransformationAnnotation ta ";

		queryStr += " join Mapping2Filter f1 on (f1.mapping2_id=mol1 and f1.filter_id=:fId)  join Mapping2Filter f2 on (f2.mapping2_id=mol2 and f2.filter_id=:fId)";

		if (mmpQuery.transformationId != null)
			queryStr += " and transformation_id=" + mmpQuery.transformationId;

		if (mmpQuery.significantTransformationsOnly)
			queryStr += " and ta.property_id=:propertyId";

		if (mmpQuery.similarity != null)
			if (mmpQuery.similarity > 0)
				queryStr += " and similarity > :similarity";
			else
				queryStr += " and similarity < -:similarity";

		ThreadScope.setStatus("Querying for MMPs for query: " + mmpQuery);

		Query q = Globals.session().createSQLQuery(queryStr)
				.addScalar("transformation_id", LongType.INSTANCE)
				.addScalar("mmp_id", LongType.INSTANCE)
				.addScalar("mol1", IntegerType.INSTANCE)
				.addScalar("mol2", IntegerType.INSTANCE)
				.setParameter("fId", mmpQuery.filterId)
				//.setParameter("propertyId", mmpQuery.propertyId)
				.setTimeout(500);

		if (mmpQuery.similarity != null)
			q.setParameter("similarity", mmpQuery.similarity);

		List<Object[]> rows =  q.list();

		if (mmpQuery.maxPairs != null && rows.size() > mmpQuery.maxPairs)
		{
			// Too many pairs. Try to gradually limit by similarity
			logger.info("Too many pairs: " + rows.size() + ". Trying to limit by similarity");
			if (mmpQuery.similarity == null)
				mmpQuery.similarity = 50;
			else if (mmpQuery.similarity <= 60)
				mmpQuery.similarity = 60;
			else if (mmpQuery.similarity <= 70)
				mmpQuery.similarity = 70;
			else if (mmpQuery.similarity <= 75)
				mmpQuery.similarity = 75;
			else if (mmpQuery.similarity <= 80)
				mmpQuery.similarity = 80;
			else
				mmpQuery.similarity = 90;

			rows = null;
			return queryByFilter(mmpQuery);
		}

		for (Object[] row : rows)
			mmps.addPair((Long) row[0], (Long) row[1], (Integer) row[2], (Integer) row[3]);

		if (mmpQuery.requestData) 
		{
			Map<Integer, FuzzyValue> molVals = MMPStatsService.getMolValues(mmpQuery.basketId, mmpQuery.propertyId, mmps.getMoleculeIds());
			mmps.molValues.put(mmpQuery.propertyId, molVals);

			if (mmpQuery.transformationId != null && mmpQuery.propertyChangeDirection != null)
			{
				List<Double> vals =  MMPStatsService.getPairDeltas(mmps, mmpQuery.transformationId, molVals);
				if (vals.size() == mmps.pairs.size())
				{
					int i = 0;
					while (i < vals.size())
					{
						if (vals.get(i) * mmpQuery.propertyChangeDirection <= 0)
						{
							vals.remove(i);
							mmps.pairs.remove(i);
							mmps.getThinPairs(mmpQuery.transformationId).remove(i);
						} else
							i++;
					}
				}
			}
		}	



		return mmps;
	}

	public MMPSubset getPairsByProperty(MMPQuery mmpQuery) {
		String key = "p-" + mmpQuery.propertyId;
		if (mmpQuery.significantTransformationsOnly)
			key += "-sig"; 
		if (mmpQuery.affectedPairsOnly)
			key += "-aff";
		if (mmpQuery.transformationId != null && !queryCache.containsKey(key))
			key += "-t-" + mmpQuery.transformationId;
		if (mmpQuery.similarity != null)
			key += "-sim-" + mmpQuery.similarity;
		if (mmpQuery.propertyChangeDirection != null)
			key += "-pcd-"+mmpQuery.propertyChangeDirection;

		if (queryCache.containsKey(key))
			return queryCache.get(key);

		mmpQuery.filterId = createMoleculeFilter(mmpQuery);
		MMPSubset mmps = queryByFilter(mmpQuery);

		Mapping2Filter.clear(mmpQuery.filterId);

		addToCache(key, mmps);

		return mmps;
	}

	public MMPSubset getPairsByMolecules(MMPQuery query) 
	{
		ThreadScope.setStatus("Querying matched pairs for a list of " + query.molIDs.size() + " molecules");

		query.filterId = Mapping2Filter.generateFilterID();
		Mapping2Filter.addCompoundToFilter(query.filterId, query.molIDs);

		MMPSubset mmps = queryByFilter(query);
		Mapping2Filter.clear(query.filterId);

		return mmps;
	}

	public MMPSubset getPairsFromPrediction(MMPQuery query) {
		String key = "applier-" + query.applier;

		if (query.insideAD)
			key += "-insideAD";

		if (query.similarity != null && query.similarity > 0)
			key += "-sim-" + query.similarity;

		if (queryCache.containsKey(key))
			return queryCache.get(key);

		ModelMapping mm = query.applier.modelTasks.get(0).model.modelMappings.get(0);
		Property property = mm.property;

		ThreadScope.setStatus("Retrieving prediction results");
		// Retrieve the prediction results for each model
		List<PredictionResults> results = new ArrayList<PredictionResults>();
		for (ModelApplierTaskProcessor mTask : query.applier.modelTasks)
		{
			for (int i = 0; i < mTask.model.modelMappings.size(); i++)
			{
				PredictionResults result = new PredictionResults();
				result.removeErrors = false;
				ApplicabilityDomain ad;
				try
				{
					ad = mTask.getApplicabilityDomain(null, 0);
				} catch (Exception e)
				{
					e.printStackTrace();
					ad = null;
				}
				result.setPredictions(mTask.wndResult.ports.get(0), ad, mTask.model.modelMappings.get(i), i);
				results.add(result);
			}
		}

		// Group predictions by molecules
		query.applier.compoundsProvider.loadStructures = false;
		Map<Integer, PropertyPrediction[]> predictionsByMols = groupPredictions(results, query.applier.compoundsProvider.getBasket());

		Map<Integer, FuzzyValue> mp2Vals = new HashMap<Integer, FuzzyValue>();
		for (Integer molId : predictionsByMols.keySet())
		{
			PropertyPrediction pp = predictionsByMols.get(molId)[0];
			if (pp != null)
			{
				if (!query.insideAD || Boolean.TRUE.equals(pp.getInsideAD()))
					mp2Vals.put(molId, new FuzzyValue(pp.getValue(), pp.getAccuracy()));
			}
		}

		if (property.isQualitative())
		{
			ThreadScope.setStatus("Resolving labels by prediction values");
			for (Integer mol : mp2Vals.keySet())
			{
				FuzzyValue val = mp2Vals.get(mol);
				val.value = mm.model.attachment.getObject().getOptionFromPrediction(val.value, mm.property).id.doubleValue();
			}
		}


		query.molIDs = predictionsByMols.keySet();

		MMPSubset subset = getPairsByMolecules(query);
		subset.predictions = predictionsByMols;
		subset.molValues.put(property.id, mp2Vals);
		addToCache(key, subset);

		return subset;
	}

	public MMPSubset getPairsFromPredictedBasket(MMPQuery mmpQuery) {

		String key = "model-" + mmpQuery.model.id + "-set-" + mmpQuery.basketId;

		if (queryCache.containsKey(key))
			return queryCache.get(key);

		ThreadScope.setStatus("Retrieving prediction results");
		// Retrieve the prediction results for each model
		List<PredictionResults> results = new ArrayList<PredictionResults>();
		int i = 0;
		for (ModelMapping mm : mmpQuery.model.modelMappings)
		{
			ModelStatistics ms = (ModelStatistics) mm.statisticsOriginal.getObject();
			PredictionResults result = new PredictionResults();
			result.removeErrors = false;
			result.setPredictions(ms.getSetStatisticsByBasket(mmpQuery.basketId), null, mm, i++);
			results.add(result);
		}

		MMPSubset subset = getPairsByBasket(mmpQuery);
		subset.predictions = groupPredictions(results, Basket.getById(mmpQuery.basketId));

		addToCache(key, subset);

		return subset;
	}

	public MMPSubset getCachedSubset(String key) {
		return queryCache.get(key);
	}




	public List<MMPTransformation> getTransformations(TransformationFilter filter)
	{
		return getTransformationsPage(filter);
	}

	private Criteria getIntermediateCriteria(TransformationFilter filter, boolean inverted) {
		Criteria c = getCriteriaFor(filter, inverted)
				.setProjection(Projections.projectionList().add(Projections.distinct(Projections.id())).add(Projections.property("pairsCount")));
		return c;
	}


	public List<MMPTransformation> getTransformationsPage(TransformationFilter filter)
	{
		List<Long> ids;
		List<Object[]> mmptIdsInverse = null;

		if (!filter.getPropertyIds().isEmpty())
		{
			List<Object[]> mmptIdsStraight = getIntermediateCriteria(filter, false).list();
			mmptIdsInverse = getIntermediateCriteria(filter, true).list();
			ids = mergeIds(mmptIdsStraight, mmptIdsInverse);
		} else
			ids = getCriteriaFor(filter, false).setProjection(Projections.projectionList().add(Projections.distinct(Projections.id()))).addOrder(Order.desc("pairsCount")).list();

		if (filter.pageInfo != null)
		{
			filter.pageInfo.size = ids.size();

			if (filter.pageInfo.pageNum > 0 && filter.pageInfo.pageSize > 0)
				ids = ids.subList((filter.pageInfo.pageNum - 1) * filter.pageInfo.pageSize, Math.min(ids.size(), filter.pageInfo.pageNum * filter.pageInfo.pageSize));
		}

		List<MMPTransformation> list = fetchByIdList(ids);

		if (mmptIdsInverse != null)
		{
			Set<Long> inverseSet = new HashSet<Long>();
			for (Object[] inv : mmptIdsInverse)
				inverseSet.add((Long)inv[0]);

			for (MMPTransformation mmpTransformation : list)
				if (inverseSet.contains(mmpTransformation.id))
					mmpTransformation.inversed = true;
		}

		return list;
	}

	private List<Long> mergeIds(List<Object[]> a, List<Object[]> b) //Merge lists maintaining sorting by count
	{
		Set<Long> s = new HashSet<Long>();
		List<Long> r = new ArrayList<Long>();
		int ai = 0, bi = 0;
		while (ai < a.size() || bi < b.size())
		{
			if (ai == a.size())
			{
				Long id = (Long)b.get(bi)[0];
				if (!s.contains(id))
				{
					r.add(id);
					s.add(id);
				}
				bi++;			
			} else
				if (bi == b.size())	
				{
					Long id = (Long)a.get(ai)[0];
					if (!s.contains(id))
					{
						r.add(id);
						s.add(id);
					}
					ai++;				
				} else
					if ((Long)a.get(ai)[1] >= (Long)b.get(bi)[1])
					{
						Long id = (Long)a.get(ai)[0];
						if (!s.contains(id))
						{
							r.add(id);
							s.add(id);
						}
						ai++;
					} else
					{
						Long id = (Long)b.get(bi)[0];
						if (!s.contains(id))
						{
							r.add(id);
							s.add(id);
						}
						bi++;
					}
		}
		return r;
	}

	private List<MMPTransformation> fetchByIdList(List<Long> mmptIds)
	{
		Map<Long, MMPTransformation> map = new HashMap<Long, MMPTransformation>();
		List<Long> ids = new ArrayList<Long>();
		ids.addAll(mmptIds);
		while (ids.size() > 0)
		{
			List<Long> subIds = ids.subList(0, Math.min(1000, ids.size()));
			List<MMPTransformation> l = Globals.session().createCriteria(MMPTransformation.class).add(Restrictions.in("id", subIds)).list();
			for (MMPTransformation mmpTransformation : l)
				map.put(mmpTransformation.id, mmpTransformation);
			subIds.clear();
		}
		List<MMPTransformation> result = new ArrayList<MMPTransformation>();
		for (Long id : mmptIds)
			result.add(map.get(id));
		return result;
	}

	private void addEffectRestriction(String alias, boolean inverted, boolean[] effects, Property p, Map<Long, Junction> criteriaMap, TransformationFilter filter) {
		Disjunction d = Restrictions.disjunction();
		int[] nums = inverted ? new int[]{2, 0} : new int[]{0, 2};
		if (effects[nums[0]])
			if (p.isQualitative())
				d.add(Restrictions.ltProperty(alias + ".subsetStats.nNP", alias + ".subsetStats.nPN"));
			else
				d.add(Restrictions.lt(alias + ".subsetStats.deltaMean", 0.0));
		if (effects[1])
			d.add(Restrictions.isNull(alias + ".id"));
		if (effects[nums[1]])
			if (p.isQualitative())
				d.add(Restrictions.geProperty(alias + ".subsetStats.nNP", alias + ".subsetStats.nPN"));
			else
				d.add(Restrictions.ge(alias + ".subsetStats.deltaMean", 0.0));

		Junction j = criteriaMap.get(p.id);
		if (j == null)
			criteriaMap.put(p.id, j = (filter.splitSets ? Restrictions.conjunction() : Restrictions.disjunction()));
		j.add(d);
	}


	private Criteria getCriteriaFor(TransformationFilter filter, boolean inverted) 
	{
		Criteria c = Globals.session().createCriteria(MMPTransformation.class)
				.addOrder(Order.desc("pairsCount"));

		if (filter.minCount != null)
			c.add(Restrictions.ge("pairsCount", filter.minCount));

		if (filter.exactCount != null)
			c.add(Restrictions.eq("pairsCount", filter.exactCount));

		if (!filter.propertyEffect.isEmpty())
			Globals.setMarshallingOption(MarshallingOption.PAIR_TRANSFORMATION);

		if (filter.ids != null)
			c.add(Restrictions.in("id", filter.ids));

		if (filter.moleculeId != null)
		{
			// Filter the transformations relevant for particular molecule
			List<Long> fragments = Globals.session().createQuery("select fragmentId from MMPIndex where mapping2Id=:mp2Id")
					.setParameter("mp2Id", filter.moleculeId)
					.list();


			c.add(Restrictions.in(inverted ? "frag2Id" : "frag1Id", fragments));
		}

		if (filter.fragId != null)
			c.add(Restrictions.or(Restrictions.eq("frag1Id", filter.fragId), Restrictions.eq("frag2Id", filter.fragId)));

		int i = 0;
		Map<Long, Junction> criteriaMap = new HashMap<Long, Junction>();
		for (Long propId : filter.getPropertyIds())
		{
			Property p = Property.getById(propId);
			if (!filter.splitSets)
			{
				Disjunction setsOr = Restrictions.disjunction();
				String alias =  "a" + propId;
				for (Long setId : filter.getSetsByProperty(propId))
					if (setId != null && setId > 0)
						setsOr.add(Restrictions.eq(alias + ".annotationSet", (MMPAnnotationSet)Globals.session().get(MMPAnnotationSet.class, setId)));
					else
						setsOr.add(Restrictions.isNull(alias + ".annotationSet"));

				c.createAlias("annotations" + (++i), alias, JoinType.LEFT_OUTER_JOIN, Restrictions.and(Restrictions.eq(alias + ".property", p), setsOr));

				boolean[] effects = filter.getEffects(propId);
				addEffectRestriction(alias, inverted, effects, p, criteriaMap, filter);
			}
			else
			{
				for (Long setId : filter.getSetsByProperty(propId))
				{
					MMPAnnotationSet set = (MMPAnnotationSet)Globals.session().get(MMPAnnotationSet.class, setId);
					Map<Long, boolean[]> effectMap = filter.propertyEffect.get(setId);
					if (i > 6)
						throw new UserFriendlyException("Too many simultaneous criteria. This is not yet supported, but will be implemented soon.");

					String alias =  "a" + propId + "s" + setId;
					boolean[] effects = effectMap.get(propId);

					Criterion joinRestrictions;
					if (setId != null && setId > 0)
						joinRestrictions = Restrictions.and(Restrictions.eq(alias + ".property", p), Restrictions.eq(alias + ".annotationSet", set));
					else
						joinRestrictions = Restrictions.and(Restrictions.eq(alias + ".property", p), Restrictions.isNull(alias + ".annotationSet"));

					c.createAlias("annotations" + (++i), alias, JoinType.LEFT_OUTER_JOIN, joinRestrictions);

					addEffectRestriction(alias, inverted, effects, p, criteriaMap, filter);
				}
			}
		}

		for (Junction j : criteriaMap.values())
			c.add(j);

		return c;
	}

	private Map<Integer, PropertyPrediction[]> groupPredictions(List<PredictionResults> results, Basket basket) {
		ThreadScope.setStatus("Grouping predictions by molecules");
		// Group predictions by molecules
		Map<Integer, PropertyPrediction[]> predictionsByMols = new HashMap<Integer, PropertyPrediction[]>();
		List<BasketEntry> entries = basket.entries;

		for (int k = 0; k < results.size(); k++)
			if (results.get(k).predictions.size() != entries.size())
				throw new UserFriendlyException("Incompatible sizes of the basket and predictions count ("+results.get(k).predictions.size()+" <> "+entries.size()+")");

		for (int i = 0; i < entries.size(); i++)
		{
			PropertyPrediction[] predictions = new PropertyPrediction[results.size()];
			for (int k = 0; k < predictions.length; k++)
				predictions[k] = results.get(k).predictions.get(i); // FIXME: Handle skipped errors

			predictionsByMols.put(entries.get(i).ep.molecule.mapping2.id, predictions);
		}

		return predictionsByMols;
	}

	private static void addToCache(String key, MMPSubset subset) {
		queryCache.put(key, subset);
		subset.id = key;

		logger.info(String.format("Identified %d MMPs for key %s", subset.pairsCount, key));

		// Some empirical cleanups of the cache
		if (MemoryUtils.getCurrentMemoryUsedFraction() > 0.8)
		{
			logger.warn("Clearing the cache because of high memory consumption");
			queryCache.clear();
		}
	}

	private static final Logger logger = LogManager.getLogger(MMPQueryService.class);

	public static void startCacheUpdate(){
		if(count == null){
			count = 10;
			new MMPQueryService().start(); // to start only one thread!
		}
		count = 10;
	}
	
	@Override
	public void run() {

		try{
			while(count != null){
				sleep(600*1000);
				if(count-- > 0)
					clearQueryCache();
			}
		}catch(InterruptedException e){}

	}

}
