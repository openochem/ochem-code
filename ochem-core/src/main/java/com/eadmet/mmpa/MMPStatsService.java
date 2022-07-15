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
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Query;
import org.hibernate.criterion.Restrictions;
import org.hibernate.type.DoubleType;
import org.hibernate.type.IntegerType;
import org.hibernate.type.LongType;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.entities.Model;
import qspr.entities.ModelAttachment;
import qspr.entities.ModelMapping;
import qspr.entities.Property;
import qspr.entities.Unit;
import qspr.modelling.ModelStatistics;
import qspr.modelling.PointStatistics;
import qspr.modelling.SetStatistics;
import qspr.util.SmartLogger;
import qspr.util.unitconversion.UnitConversion;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.mmpa.domain.FuzzyValue;
import com.eadmet.mmpa.domain.MMPSubset;
import com.eadmet.mmpa.domain.MMPSubsetStatistics;
import com.eadmet.mmpa.domain.MMPTransformation;
import com.eadmet.mmpa.domain.ThinPair;
import com.eadmet.utils.MemCache;

/**
 * A service for calculations various statistical parameters of MMPs
 * 
 * @author midnighter
 *
 */
@SuppressWarnings("unchecked")
public class MMPStatsService
{
	public static List<MMPTransformation> getTransformationsStats(MMPStatsRequest request) {

		MMPQueryService service = new MMPQueryService();
		Map<Integer, FuzzyValue> mp2valReal;
		MMPSubset subset;

		MMPQuery query = new MMPQuery();
		query.setSimilarity(request.similarity);
		query.maxPairs = request.maxPairs;

		if (request.subset != null) {
			subset = request.subset;
			mp2valReal = request.subset.molValues.get(request.subset.molValues.keySet().iterator().next());
		} else if (request.applier != null) {
			request.property = request.applier.modelTasks.get(0).model.modelMappings.get(0).property;
			request.property = Property.getById(request.property.id);

			query.applier = request.applier;
			query.insideAD = request.insideAD;
			subset = service.getPairsFromPrediction(query);
			mp2valReal = subset.molValues.values().iterator().next();
		}
		else if (request.basket != null) 
		{
			subset = service.getPairsByBasket(query.setBasket(request.basket.id).setProperty(request.property.id));
			mp2valReal = getMolValues(request.basket.id, request.property.id, null);
		}
		else if (request.property != null && request.model == null) 
		{
			subset = service.getPairsByProperty(query.setProperty(request.property.id));
			mp2valReal = getMolValues(null, request.property.id, null);
		}
		else
		{
			if(request.property == null)
				request.property = request.model.modelMappings.get(0).property;
			query.propertyId = request.property.id;
			subset = service.getPairsByBasket(query.setBasket(request.model.trainingSet.id));
			mp2valReal = getMolValues(request.model, request.usePredictedValues, request.property);
		}

		if (request.property != null || request.model != null)
		{
			Property p = request.property != null ? request.property : request.model.modelMappings.get(0).property;
			// Cache the values to nicely display them in UI
			subset.molValues.put(p.id, mp2valReal);
		}

		logger.info("Loaded experimental values for " + mp2valReal.size() + " molecules");

		Map<Long, Long> labelsMapping = null;
		boolean classification = request.property != null && request.property.isQualitative();
		if (classification)
			labelsMapping = Model.getClassificationRemapping(request.property);
		request.subset = subset;

		List<MMPTransformation> chosenTransformations = new ArrayList<MMPTransformation>();

		List<Long> transformationIds = new ArrayList<Long>();
		for (Long trId : subset.getTransformationIds())
			if (subset.getThinPairs(trId).size() >= request.minPairs)
				transformationIds.add(trId);

		if (transformationIds.isEmpty())
			return chosenTransformations;

		// Preloading transformations for efficiency
		Globals.session().createCriteria(MMPTransformation.class).add(Restrictions.in("id", transformationIds)).list();

		Map<Long, Integer> atomCounts = null;
		if (request.maxAtoms != null)
			atomCounts = getTransformationSizes(transformationIds);

		int bootstrapReplicas = request.bootstrapStatistics ? 100 : 1;
		MersenneTwister random = new MersenneTwister();

		Map<Long, IntermediateStats> cachedStats = statsCache.get(request.getMinorKey());
		if (cachedStats == null)
			statsCache.put(request.getMinorKey(), cachedStats = new HashMap<Long, IntermediateStats>());
		// Calculate stats for every prospective transformation
		ThreadScope.setStatus("Calculating stats for the transformations", logger);
		SmartLogger sLogger = new SmartLogger(MMPStatsService.logger, 1000);
		if (ThreadScope.get().operation != null)
			sLogger.setStatusTracker(ThreadScope.get().operation.statusTracker);
		int cnt = 0;

		String []fragments = null;

		if(request.transformationPattern != null){
			fragments = request.transformationPattern.split(",");
		}

		for (Long trId : transformationIds) {

			cnt++;
			sLogger.log("" + cnt + " transformations out of " + transformationIds.size());
			IntermediateStats is = cachedStats.get(trId);
			if (is == null)
			{
				is = new IntermediateStats();

				for (int r = 0; r < bootstrapReplicas; r++) 
				{
					List<ThinPair> pairs = subset.getThinPairs(trId);

					for (ThinPair pair : pairs)
					{
						FuzzyValue[] fv = new FuzzyValue[]{mp2valReal.get(pair.mol1Id), mp2valReal.get(pair.mol2Id)};
						if (fv[0] == null || fv[1] == null)
							continue;

						double delta = 0;
						if (classification)
						{
							long[] val = new long[]{labelsMapping.get(Double.valueOf(fv[0].value).longValue()), labelsMapping.get(Double.valueOf(fv[1].value).longValue())};

							// Perturbate the prediction
							if (request.bootstrapStatistics)
								for (int k = 0; k < 2; k++)
									if (fv[k].accuracy != null && random.nextInt(1024) >= 1024 * fv[k].accuracy)
										val[k] = - val[k];

							delta = val[1] - val[0];
							if (val[1] == val[0])
								if (val[1] == -1)
									is.nEqNN++;
								else
									is.nEqPP++;
						}
						else
						{
							double[] val = new double[]{mp2valReal.get(pair.mol1Id).value, mp2valReal.get(pair.mol2Id).value};

							// Perturbate the prediction
							if (request.bootstrapStatistics)
								for (int k = 0; k < 2; k++)
									if (fv[k].accuracy != null)
										val[k] = val[k] + random.nextGaussian() * fv[k].accuracy;
							delta = val[1] - val[0];
						}

						if (delta > 0)
							is.nPos++;
						else if (delta < 0)
							is.nNeg++;


						is.sum += delta;
						is.sum2 += delta * delta;
						is.count++;
					}
				}
				cachedStats.put(trId, is);
			}

			if (is.count >= request.minPairs * bootstrapReplicas)
			{
				double pValue = 0;
				if (classification)
					pValue = is.getPValueClassification(bootstrapReplicas);
				else
					pValue = new BinomialDistribution((is.nPos + is.nNeg) / bootstrapReplicas, 0.5).cumulativeProbability(Math.min(is.nPos / bootstrapReplicas, is.nNeg / bootstrapReplicas));
				if (pValue <= request.pValue)
				{

					if (request.maxAtoms != null)
						if (atomCounts.get(trId) > request.maxAtoms)
							continue;

					MMPTransformation transformation = MMPTransformation.get(trId);

					if(transformation == null){
						System.out.println ("transformation == null " + trId + " " +is);
						return chosenTransformations;
					}

					transformation.statistics = new MMPSubsetStatistics(is.sum / is.count, is.count == 0 ? 0: 
						Math.sqrt((is.sum2 - is.sum*is.sum / is.count) / (is.count - 1)), is.count / bootstrapReplicas);
					transformation.statistics.pValue = pValue;
					if (classification)
					{
						transformation.statistics.nNN = Math.round(1.0f * is.nEqNN / bootstrapReplicas);
						transformation.statistics.nNP = Math.round(1.0f * is.nPos / bootstrapReplicas);
						transformation.statistics.nPN = Math.round(1.0f * is.nNeg / bootstrapReplicas);
						transformation.statistics.nPP = Math.round(1.0f * is.nEqPP / bootstrapReplicas);
					}

					if (!classification)
						if (request.meanStdFactor != null && Math.abs(transformation.statistics.deltaMean) < request.meanStdFactor * Math.abs(transformation.statistics.deltaStd))
							continue;

					// Simple filtering by substructure match

					if(request.transformationPattern != null){
						boolean found = false;
						for(String fragment:fragments)
							if(transformation.getSmirks().contains(fragment)){
								found = true;
								break;
							}
						if(!found)continue;
					}

					chosenTransformations.add(transformation);
				}
			}
		}
		logger.info("Statistics calculated");


		holmBolferonni(chosenTransformations, request.pValue);

		Collections.sort(chosenTransformations, new Comparator<MMPTransformation>()
		{
			@Override
			public int compare(MMPTransformation arg0, MMPTransformation arg1)
			{
				return new Double(arg0.statistics.deltaMean).compareTo(new Double(arg1.statistics.deltaMean));
			}
		});

		statsCache.cleanup(5);

		return chosenTransformations;
	}

	/**
	 * Correction of false positives according to Holm-Bolferonni method (see Wikipedia).
	 * This method has been suggested by one of the publication reviewers.
	 */
	private static void holmBolferonni(List<MMPTransformation> transformations, double significance) {
		if (significance >= 1.0)
			return;

		Collections.sort(transformations, new Comparator<MMPTransformation>()
		{
			@Override
			public int compare(MMPTransformation arg0, MMPTransformation arg1)
			{
				return new Double(arg0.statistics.pValue).compareTo(new Double(arg1.statistics.pValue));
			}
		});

		int i = 0;
		Iterator<MMPTransformation> iTrans = transformations.iterator();
		while (iTrans.hasNext())
		{
			if (iTrans.next().statistics.pValue > significance / (transformations.size() - i))
			{
				logger.info( "" + (transformations.size() - i) + " out of " + transformations.size() + " transformations as insignificant according to Holm-Bolferonni method");
				iTrans.remove();
				break;
			}
			i++;
		}

		while (iTrans.hasNext())
		{
			iTrans.next();
			iTrans.remove();
		}
	}

	/**
	 * Get the property values grouped by molecules
	 */
	public static Map<Integer, FuzzyValue> getMolValues(Model model, boolean predicted, Property property) {
		logger.info("Grouping " + (predicted ? "predicted" : "real") + " values by molecules for model " + model.name);

		ModelStatistics ms = null;
		if(property!=null)
			for(ModelMapping m: model.modelMappings)
				if(property == m.property)
					ms = (ModelStatistics)m.statisticsOriginal.getObject();
		if(ms == null){
			ms = (ModelStatistics) model.modelMappings.get(0).statisticsOriginal.getObject();
			property = model.modelMappings.get(0).property;
		}
		SetStatistics ss = ms.getSetStatisticsByBasket(model.trainingSet.id);
		Map<Long, Integer> ep2mp2 = new HashMap<Long, Integer>();
		Map<Integer, FuzzyValue> mp2val = new HashMap<Integer, FuzzyValue>();
		if(ss == null || ss.getEpIds() == null || ss.getEpIds().size() == 0) return mp2val;
		List<Object[]> rows = Globals.session().createQuery("select ep.id, mp2.id from ExperimentalProperty ep join ep.molecule m join m.mapping2 mp2 where ep.id in (:ids)")
				.setParameterList("ids", ss.getEpIds())
				.list();

		ModelAttachment attachment = model.attachment.getObject();
		for (Object[] row : rows)
			ep2mp2.put((Long) row[0], (Integer) row[1]);
		for (PointStatistics ps : ss.points)
		{
			if (predicted && ps.error != null)
				continue;
			if (!property.isQualitative())
				mp2val.put(ep2mp2.get(ps.id), new FuzzyValue(predicted ? ps.predicted : ps.real));
			else
				mp2val.put(ep2mp2.get(ps.id), new FuzzyValue(attachment.getOptionFromPrediction(predicted ? ps.predicted : ps.real, property).id.doubleValue()));
		}

		return mp2val;
	}


	public static Map<Integer, FuzzyValue> getMolValues(Long basketId, long propertyId, Collection<Integer> moleculeIds) {

		Property property = Property.getById(propertyId);
		boolean classification = property.isQualitative();

		String valueColumn = classification ? "poption_id" : "value";

		String queryStr = "select mapping2_id, " + valueColumn + " val, unit_id, mol_weight from ExperimentalProperty join Molecule using (molecule_id)";
		if (basketId != null)
			queryStr += " join BasketEntry using (exp_property_id)";
		queryStr += " where property_id=:propertyId";
		if (basketId != null)
			queryStr += " and basket_id=:basketId";
		if (basketId == null)
			queryStr += " and rights=" + Globals.RIGHTS_FREELY_AVAILABLE;
		if (basketId == null)
			queryStr += " and exp_property_id=first_entry";
		if (moleculeIds != null)
			if (moleculeIds.size() > 0)
				queryStr += " and mapping2_id in (:molIds)";
			else
				queryStr += " and (1 = 0)";

		queryStr += " group by mapping2_id";

		Query query = Globals.session().createSQLQuery(queryStr)
				.addScalar("mapping2_id", IntegerType.INSTANCE)
				.addScalar("val", DoubleType.INSTANCE)
				.addScalar("unit_id", LongType.INSTANCE)
				.addScalar("mol_weight", DoubleType.INSTANCE)
				.setParameter("propertyId", propertyId);

		if (moleculeIds != null)
			if (moleculeIds.size() > 0)
				query.setParameterList("molIds", moleculeIds);

		if (basketId != null)
			query.setParameter("basketId", basketId);

		List<Object[]> rows = query.list();

		Map<Integer, FuzzyValue> mp2Val = new HashMap<Integer, FuzzyValue>();

		logger.info("Unit conversion and grouping of values by molecules...");

		for (Object[] row : rows)
		{
			Double value = (Double) row[1];
			if (!classification)
			{
				Long unitId = (Long) row[2];
				if (!unitId.equals(property.defaultUnit.id))
					try
				{
						value = UnitConversion.convert(value, Unit.getById(unitId), property.defaultUnit, (Double) row[3]);
				} catch (UserFriendlyException e) {
					logger.warn(e);
				}
			}
			mp2Val.put((Integer) row[0], new FuzzyValue(value));
		}

		logger.info("Got values for " + mp2Val.size() + " molecules");


		return mp2Val;
	}

	public static List<Double> getPairDeltas(long propertyId, long transformationId) {
		MMPQuery query = new MMPQuery();
		query.propertyId = propertyId;
		query.transformationId = transformationId;
		MMPSubset subset = new MMPQueryService().getPairsByProperty(query);
		Map<Integer, FuzzyValue> molVals = getMolValues(null, propertyId, subset.getMoleculeIds());

		return getPairDeltas(subset, transformationId, molVals);
	}

	public static List<Double> getPairDeltas(String subsetKey, long transformationId) {
		MMPSubset subset = new MMPQueryService().getCachedSubset(subsetKey);
		return getPairDeltas(subset, transformationId, subset.molValues.get(subset.molValues.keySet().iterator().next()));
	}

	public static List<Double> getPairDeltas(MMPSubset subset, long transformationId, Map<Integer, FuzzyValue> molVals) {
		List<Double> deltas = new ArrayList<Double>();
		List<ThinPair> pairs = subset.getThinPairs(transformationId);
		if (pairs == null)
			return deltas;
		for (int i = 0; i < pairs.size(); i++)
		{
			ThinPair pair = pairs.get(i);
			if (molVals.containsKey(pair.mol1Id) && molVals.containsKey(pair.mol2Id))
				deltas.add(molVals.get(pair.mol2Id).value - molVals.get(pair.mol1Id).value);
		}

		return deltas;	
	}

	private static Map<Long, Integer> getTransformationSizes(Collection<Long> trIds) {
		Map<Long, Integer> map = new HashMap<Long, Integer>();
		List<Object[]> rows = Globals.session().createSQLQuery("select transformation_id tid, f1.size s1, f2.size s2 from MMPTransformation join MMPFragment f1 on (f1.fragment_id=fragment1_id) join MMPFragment f2 on (f2.fragment_id=fragment2_id) where (transformation_id in (:ids))")
				.addScalar("tid", LongType.INSTANCE)
				.addScalar("s1", IntegerType.INSTANCE)
				.addScalar("s2", IntegerType.INSTANCE)
				.setParameterList("ids", trIds)
				.list();

		for (Object[] row : rows)
			map.put((Long)row[0], Math.max((Integer) row[1], (Integer) row[2]));

		return map;
	}

	public static void main(String[] args)
	{
		IntermediateStats is = new IntermediateStats();
		is.nNeg = 0;
		is.nEqNN = 0;
		is.nEqPP = 1;

		System.out.println(is.getPValueClassification(1));
	}

	private static final Logger logger = LogManager.getLogger(MMPStatsService.class);
	private static MemCache<String, Map<Long, IntermediateStats>> statsCache = new MemCache<String, Map<Long,IntermediateStats>>();
}

class IntermediateStats {
	double sum = 0, sum2 = 0;
	int nPos = 0, nNeg = 0, nEqPP = 0, nEqNN = 0;
	int count = 0;


	public String toString(){
		return "sum = " + sum + " sum2 = " +sum2 + " count = " + count;
	}

	/**
	 * A sophisticated calculation of pValue for classification properties.
	 * Our "null hypothesis" is a binomial distribution of positives according to the ration for molecules BEFORE (or after) the transformation
	 */
	public double getPValueClassification(int bootstrapReplicas) {

		// Transformation effect counts (00, 01, 10, and 11)
		int[][] effCount = new int[][]{{nEqNN, nPos}, {nNeg, nEqPP}};

		int totals = effCount[0][0] + effCount[0][1] + effCount[1][0] + effCount[1][1];
		int[] positives = new int[]{effCount[1][0] + effCount[1][1], effCount[0][1] + effCount[1][1]};

		// Special singular case
		if (positives[0] == 0 && positives[1] == 0)
			return 0.5;

		if (positives[0] == positives[1])
			return 0.5;

		// Choose standard (null-hypothesis) distribution among BEFORE and AFTER the transformation. 
		// Namely, choose the one having the positives ratio closer to 0.5
		int standard = Math.abs(2 * positives[0] - totals) < Math.abs(2 * positives[1] - totals) ? 0 : 1;

		if (positives[standard] == totals)
			standard = 1 - standard;

		if (positives[standard] == 0)
		{
			positives[standard] = 1;
			totals++;
		}

		if (positives[1 - standard] == totals)
			positives[1 - standard]--;

		// Divide by bootstrap replicas to avoid "significance over-estimation" due to bootstrapping
		double pValue = new BinomialDistribution(totals / bootstrapReplicas, 1.0 * positives[standard] / totals).cumulativeProbability(positives[1 - standard] / bootstrapReplicas);
		if (pValue > 0.5)
			pValue = 1 - pValue;
		return pValue;
	}
}
