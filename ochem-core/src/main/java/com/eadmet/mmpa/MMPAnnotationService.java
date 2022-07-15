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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.sql.JoinType;
import org.hibernate.type.LongType;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.Property;
import qspr.util.WrapperThread;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.mmpa.domain.MMPAnnotationSet;
import com.eadmet.mmpa.domain.MMPTransformation;
import com.eadmet.mmpa.domain.MMPTransformationAnnotation;

@SuppressWarnings("unchecked")
public class MMPAnnotationService
{
	public static AnnotationCountsCache countsCache = new AnnotationCountsCache();

	public static void indexMolecules(Basket basket){
		int mol = 0;
		try{
			mol = MMPIndexingService.getInstance().reSubmitIndexingTasks(basket);
		}catch(Exception e){
		}
		throw new UserFriendlyException( mol == 0? "All molecules have been alraedy indexed and no new molecules are submitted.":"Molecules N= "+ mol+" are submited. Results will be ready a couple of hours (depending on the server load). ");
	}

	public static void annotateRequestedTransformations(MMPStatsRequest request) 
	{
		annotateRequestedTransformations(request, null); 
	}

	private static Criteria getCriteria(MMPAnnotationSet aset)
	{
		Criteria c = Globals.session().createCriteria(MMPTransformationAnnotation.class);
		c.createAlias("annotationSet", "aset", JoinType.LEFT_OUTER_JOIN);

		List<Long> ids = aset.getIdList();
		Disjunction d = Restrictions.disjunction();
		for (Long id : ids) 
		{
			if (id > 0)
				d.add(Restrictions.eq("aset.id", id));
			else
				d.add(Restrictions.isNull("aset.id"));
		}
		c.add(d);
		return c;
	}

	public static void deleteInvalidAnnotations() {
		List<Long> annosToDelete = Globals.session().createSQLQuery("select tp_id from MMPTransformationAnnotation left join MMPTransformation t using (transformation_id) where t.transformation_id is null").addScalar("tp_id", LongType.INSTANCE).list();
		if (!annosToDelete.isEmpty())
			Globals.session().createSQLQuery("delete from MMPTransformationAnnotation where tp_id in (:ids)").setParameterList("ids", annosToDelete).executeUpdate();
	}

	public static void calculateAnnotationCounts(MMPAnnotationSet aset) 
	{
		List<Property> properties = getCriteria(aset).setProjection(Projections.groupProperty("property")).list();
		aset.properties = properties;

		boolean allCached = true;
		for (Property property : properties) 
		{
			Globals.session().evict(property);
			if (countsCache.getCachedCounts(aset.getIdList(), property.id) != null)
				property.transformationsCount = countsCache.getCachedCounts(aset.getIdList(), property.id);
			else
				allCached = false;
		}

		if (allCached)
			return;

		Map<Property, Integer> index = new HashMap<Property, Integer>();
		for (int i = 0; i < properties.size(); i++)
			index.put(properties.get(i), i);

		int[][] counts = new int[properties.size()][2];
		List<MMPTransformationAnnotation> annotations = getCriteria(aset).list();

		Set<String> uniqueTransformations = new HashSet<String>(); //Distinct transformations within each property
		for (MMPTransformationAnnotation annotation : annotations)
		{


			boolean increasing;
			if (annotation.property.isQualitative())
				increasing = annotation.subsetStats.nNP >= annotation.subsetStats.nPN;
				else
					increasing = annotation.subsetStats.deltaMean >= 0;

					String key = annotation.property.id + "_" +annotation.transformation.id + "_" + increasing;

					if (!uniqueTransformations.contains(key))
					{
						counts[index.get(annotation.property)][increasing ? 1 : 0] ++;
						uniqueTransformations.add(key);
					}
		}

		for (int i = 0; i < properties.size(); i++)
		{
			properties.get(i).transformationsCount = counts[i];
			countsCache.putCachedCounts(aset.getIdList(), properties.get(i).id, counts[i]);
		}
	}

	public static void annotateRequestedTransformations(MMPStatsRequest request, String annotationSetName) 
	{
		List<MMPTransformation> transformations = MMPStatsService.getTransformationsStats(request);

		MMPAnnotationSet set = null;
		if (annotationSetName != null)
			set = MMPAnnotationSet.getByName(annotationSetName);

		if (set != null)
			if (Globals.myself() == null || !Globals.myself().equals(set.user))
				throw new UserFriendlyException("You cannot change annotation set of user " + set.user.login);

		int cnt = 0;
		logger.info("Deleting duplicate annotations");

		List<Long> transformationIds = new ArrayList<Long>();
		for (MMPTransformation transformation : transformations)
			transformationIds.add(transformation.id);

		if (set != null)
		{
			Globals.session().createSQLQuery("delete from MMPTransformationAnnotation where as_id=:as_id and property_id=:prop_id and transformation_id in :trans_id")
			.setLong("as_id", set.id)
			.setLong("prop_id", request.property.id)
			.setParameterList("trans_id", transformationIds)
			.executeUpdate();

			countsCache.invalidateCachedCounts(set.id, request.property.id);
		} else
		{
			Globals.session().createSQLQuery("delete from MMPTransformationAnnotation where as_id is null and property_id=:prop_id and transformation_id in :trans_id")
			.setLong("prop_id", request.property.id)
			.setParameterList("trans_id", transformationIds)
			.executeUpdate();

			countsCache.invalidateCachedCounts(0L, request.property.id);
		}

		logger.info("Starting annotating transformations");
		for (MMPTransformation transformation : transformations)
		{
			MMPTransformationAnnotation annotation = new MMPTransformationAnnotation();
			annotation.transformation = transformation;
			annotation.property = request.property;
			annotation.annotationSet = set;
			annotation.subsetStats = transformation.statistics;

			if (request.basket != null)
				annotation.basket = request.basket;

			if (request.model != null)
				annotation.model = request.model;

			if (request.applier != null)
				annotation.model = request.applier.modelTasks.get(0).model;

			Globals.session().save(annotation);
			cnt++;
			if (cnt % 1000 == 0)
			{
				logger.info("So far "+cnt+" transformations annotated");
				Globals.restartAllTransactions(true);
			}
		}
		logger.info("" + cnt + " transformations annotated for request "+request);
	}


	public static void annotateProperty(Property p) 
	{
		logger.info("Annotating transformations for property " + p.getName());
		MMPStatsRequest request = new MMPStatsRequest();
		request.minPairs = 10;
		request.pValue = 0.01;
		request.property = p;
		request.publicData = true;
		request.meanStdFactor = 0.8;
		request.primaryRecords = true;
		Globals.session().createQuery("delete from MMPTransformationAnnotation where property=:property and annotationSet is null")
		.setParameter("property", request.property)
		.executeUpdate();

		annotateRequestedTransformations(request);
	}

	public static void main(String[] args)
	{
		new WrapperThread()
		{
			@Override
			public void wrapped() throws Exception
			{
				String[] propNames = new String[]{"AMES"};
				//new String[]{"log(IGC50-1)"};
				for (String propName : propNames)
				{
					annotateProperty(Property.getByName(propName));
					Globals.restartAllTransactions(true);
				}
			}
		}.run();
	}

	private static final Logger logger = LogManager.getLogger(MMPAnnotationService.class);
}

class AnnotationCountsCache
{
	class CacheKey
	{
		@Override
		public int hashCode() 
		{
			return toString().hashCode();
		}

		@Override
		public boolean equals(Object o)
		{
			if (!(o instanceof CacheKey))
				return false;
			CacheKey key = (CacheKey) o;
			return toString().equals(key.toString());
		}

		List<Long> setIds;
		Long propertyId;

		public CacheKey(List<Long> setIds, Long propertyId)
		{
			this.setIds = setIds;
			this.propertyId = propertyId;
		}

		public boolean equalSetIds(List<Long> setIds)
		{
			Collections.sort(setIds);
			return (setIds.toString().equals(this.setIds.toString()));
		}

		public boolean containsSetId(Long setId)
		{
			return setIds.contains(setId);
		}

		public String toString()
		{
			return "{"+setIds.toString() + ", " + propertyId+"}";
		}
	}

	Map<CacheKey, int[]> cachedCounts = new HashMap<CacheKey, int[]>();

	public void clear() {
		cachedCounts.clear();
	}

	public void invalidateCachedCounts(Long setId, Long propertyId)
	{
		List<CacheKey> c = new ArrayList<CacheKey>();
		for (CacheKey key : cachedCounts.keySet())
			if (key.containsSetId(setId))
			{
				logger.info("Invalidating key = "+key);
				c.add(key);
			}
		for (CacheKey cacheKey : c)
			cachedCounts.remove(cacheKey);
	}

	public void invalidateCahcedCounts(List<Long> setIds, Long propertyId)
	{
		if (propertyId != null)
		{
			CacheKey key = new CacheKey(setIds, propertyId);
			logger.info("Invalidating key = "+key);
			cachedCounts.remove(key);
		}
		else // To invalidate for example all Default counts
		{
			List<CacheKey> c = new ArrayList<CacheKey>();
			for (CacheKey key : cachedCounts.keySet())
				if (key.equalSetIds(setIds))
				{
					logger.info("Invalidating key = "+key);
					c.add(key);
				}
			for (CacheKey cacheKey : c)
				cachedCounts.remove(cacheKey);
		}
	}

	public int[] getCachedCounts(Long setId, Long propertyId)
	{
		List<Long> setIds = new ArrayList<Long>();
		setIds.add(setId);
		return getCachedCounts(setIds, propertyId);
	}

	public void putCachedCounts(Long setId, Long propertyId, int[] counts)
	{
		List<Long> setIds = new ArrayList<Long>();
		setIds.add(setId);
		putCachedCounts(setId, propertyId, counts);
	}

	public int[] getCachedCounts(List<Long> setIds, Long propertyId)
	{
		CacheKey key = new CacheKey(setIds, propertyId);
		int[] counts = cachedCounts.get(key);
		if (counts != null)
			logger.info("Found in cache key = "+key+" value = ["+counts[0]+", "+counts[1]+"]");
		return counts;
	}

	public void putCachedCounts(List<Long> setIds, Long propertyId, int[] counts)
	{
		CacheKey key = new CacheKey(setIds, propertyId);
		logger.info("Caching key = "+key+" value = ["+counts[0]+", "+counts[1]+"]");
		cachedCounts.put(key, counts);
	}

	private static final Logger logger = LogManager.getLogger(AnnotationCountsCache.class);
}
