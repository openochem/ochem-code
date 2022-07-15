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

package qspr.frontend;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.xml.bind.annotation.XmlAnyElement;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElements;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.Query;
import org.hibernate.ScrollableResults;
import org.hibernate.criterion.Projection;
import org.hibernate.criterion.ProjectionList;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.transform.ResultTransformer;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.toxicity.RequestStatistics;
import qspr.util.CriteriaWrapper;
import qspr.workflow.utils.QSPRConstants;

@XmlRootElement(name = "list")
@SuppressWarnings({"unchecked","rawtypes"})
public class WebList 
{
	@XmlAnyElement
	public List<Object> list;

	@XmlElements
	({@XmlElement(name = "string", type=String.class)})
	public List<Object> nativeObjectsList;

	@XmlAttribute(name = "size")
	public int size;

	@XmlAttribute
	public int pageNum;

	@XmlAttribute
	public int pageSize;

	@XmlTransient 
	public String column = null;

	@XmlTransient
	private Class entityClass;

	/**
	 * We might need criteria wrapper here for some additional information like ordering and necessary aliases
	 */
	@XmlTransient
	private CriteriaWrapper criteriaWrapper;

	@XmlTransient
	public CriteriaWrapper getCriteriaWrapper()
	{
		return criteriaWrapper;
	}

	public void setCriteriaWrapper(CriteriaWrapper criteriaWrapper)
	{
		this.criteriaWrapper = criteriaWrapper;
	}

	@XmlTransient
	private Projection projection;

	public Double queryTime; // in seconds

	/**
	 * A flag that specifies that the sorting specification has been discarded for this request because of performance reasons
	 */
	public boolean sortingDiscarded;

	static Map <Integer,QueryCachedSize> hashedSize = new HashMap<Integer,QueryCachedSize>();


	/**
	 * If more than this, no sorting will be applied because of performance reasons
	 */
	public int maxAllowedItemsForSorting = 10000;

	public WebList loadFromPagedList(List list, int pageNum, int pageSize, long size)
	{
		this.list = list;
		this.size = Long.valueOf(size).intValue();
		this.pageNum = pageNum;
		this.pageSize = pageSize;
		return this;		
	}

	public WebList loadFromList(List list)
	{
		this.list = list;
		size = list.size();
		pageNum = 1;
		pageSize = list.size();
		return this;
	}

	public WebList loadFromList(List list, int pageNum, int pageSize)
	{
		this.list = list.subList((pageNum-1)*pageSize, Math.min(pageNum*pageSize,list.size()));
		size = list.size();
		this.pageNum = pageNum;
		this.pageSize = pageSize;
		return this;
	}

	/**
	 * We have IDs of all the matched elements and we want to load the respective elements 
	 * @return
	 */
	public WebList loadFromIDs(Criteria c, List ids, int pageNum, int pageSize)
	{
		if (ids.isEmpty())
			return null;
		c.add(Restrictions.in("id", ids.subList((pageNum-1)*pageSize, Math.min(pageNum*pageSize, ids.size()))));
		list = c.list();
		this.pageNum = pageNum;
		this.pageSize = pageSize;
		size = ids.size();
		return this;
	}

	public WebList loadFromSet(Set set)
	{
		this.list = new ArrayList(set);
		size = set.size();
		pageNum = 1;
		pageSize = set.size();
		return this;
	}

	public void loadFromQuery(Query query, int pageNum, int pageSize)
	{
		ScrollableResults scroll = query.scroll();
		scroll.last();
		size = scroll.getRowNumber() + 1;

		list = query.setMaxResults(pageSize).setFirstResult((pageNum-1)*pageSize).setResultTransformer(new MyEntityResultTransformer()).list();
		this.pageNum = pageNum;
		this.pageSize = pageSize;
	}

	int getSize(Criteria query){
		int hash = query.toString().hashCode();

		if(hashedSize.containsKey(hash) && !hashedSize.get(hash).isExpired(maxAllowedItemsForSorting) ){
			logger.info("Using hashed size value " + hashedSize.get(hash).size);
			return hashedSize.get(hash).size;
		}

		ProjectionList projList = Projections.projectionList();
		RequestStatistics stats = new RequestStatistics();
		stats.start();
		projList.add(Projections.count(column));
		query.setProjection(projList);
		ThreadScope.setStatus("Counting the matching results...");
		int size = ((Long) query.uniqueResult()).intValue();
		logger.info("Count query with " + size + " results executed in " + stats.current());
		hashedSize.put(hash, new QueryCachedSize(size));
		logger.info("List query executed: " + stats.current());

		logger.info("There is a new hash: "+hash + " for: " + size + " " + query.toString());
		
		return size;
	}

	public void loadFromCriteria(Criteria query, int pageNum, int pageSize)
	{
		if (column == null)
			column = "id";

		size = getSize(query);

		sortingDiscarded = size > maxAllowedItemsForSorting;

		addOrderTo(query);

		if (!column.equals("id") && ! column.equals("mapping1_id"))
			query.setProjection(Projections.groupProperty(column));
		else
			query.setProjection(null);

		if (projection == null)
			query.setResultTransformer(Criteria.ROOT_ENTITY);
		else
			query.setProjection(projection);
		ThreadScope.setStatus("Quering "+pageSize+" results out of " + size + "...");
		list = query.setMaxResults(pageSize).setFirstResult((pageNum-1)*pageSize).list();
		this.pageNum = pageNum;
		this.pageSize = pageSize;
	}

	public void loadDistinctFromCriteria(Criteria query, int pageNum, int pageSize)
	{
		// Get count of items
		ProjectionList projList = Projections.projectionList();
		projList.add(Projections.countDistinct("id"));
		query.setProjection(projList);
		long time = Calendar.getInstance().getTimeInMillis();
		ThreadScope.setStatus("Counting the matching results...");
		size = ((Long)query.uniqueResult()).intValue();
		logger.info("Count query with " + size + " results executed in " + (Calendar.getInstance().getTimeInMillis() - time) + "ms.");

		sortingDiscarded = size > maxAllowedItemsForSorting;

		// Get IDs of items on selected page
		query.setProjection(Projections.distinct(Projections.id()));
		time = Calendar.getInstance().getTimeInMillis();
		ThreadScope.setStatus("Getting the record IDs for the current page...");
		List ids = query.setMaxResults(pageSize).setFirstResult((pageNum-1)*pageSize).list();
		logger.info("List query executed: " + (Calendar.getInstance().getTimeInMillis() - time));


		this.pageNum = pageNum;
		this.pageSize = pageSize;

		if (ids.size()<=0)
		{
			list = new ArrayList();
			return;
		}

		// Get objects themselves
		Criteria helper = Globals.session().createCriteria(entityClass);
		helper.add(Restrictions.in("id", ids));

		if (projection != null)
			helper.setProjection(projection);

		// FIXME: We got a problem here. What if the order definition requires some joins?
		if (criteriaWrapper != null)
			criteriaWrapper.addAliasesTo(helper);
		addOrderTo(helper);

		helper.setResultTransformer(Criteria.ROOT_ENTITY);

		// Ok. Now post-sort the result in order, going along woth Ids list
		// What to do? Have better solution?
		// Hibernate Criteria + Pagination + OneToMany + Distinct + Sorting = Problems
		time = Calendar.getInstance().getTimeInMillis();
		ThreadScope.setStatus("Getting the record data for the current page...");
		List tempList = helper.list();

		// Remove dublicated entries (how to do it unioversally in criteria query?)
		Set set = new HashSet();
		if (!tempList.isEmpty())
		{
			set.addAll(tempList);
			tempList.clear();
			tempList.addAll(set);
		}

		list = new ArrayList();
		list.addAll(tempList);
		try {
			for (int i = 0; i < tempList.size(); i++)
			{
				Object obj = tempList.get(i);
				Long id = (Long)obj.getClass().getField("id").get(obj);
				list.set(ids.indexOf(id), obj);
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
		logger.info("Post-sorting of result finished: " + (Calendar.getInstance().getTimeInMillis() - time));
	}

	private void addOrderTo(Criteria c)
	{
		if (!sortingDiscarded && criteriaWrapper != null)
			criteriaWrapper.addOrderTo(c);
	}

	public WebList useEntity(Class classEntity)
	{
		this.entityClass = classEntity;
		return this;
	}

	public WebList useProjection(Projection projection)
	{
		this.projection = projection;
		return this;
	}

	@XmlAttribute(name = "firstResult")
	public int getFirstResult()
	{
		return (pageNum-1)*pageSize + 1;
	}

	@XmlAttribute(name = "lastResult")
	public int getLastResult()
	{
		return Math.min(size, pageNum*pageSize);
	}

	//public WebList setOrder(List<Order> order) {
	//	this.order = order;
	//	return this;
	//}

	private static transient final Logger logger = LogManager.getLogger(WebList.class);

}

class MyEntityResultTransformer implements ResultTransformer {

	private static final long serialVersionUID = 1L;

	public Object transformTuple(Object[] tuple, String[] aliases) {
		return tuple[0];
	}

	@SuppressWarnings("rawtypes")
	public List transformList(List collection) {
		return collection;
	}		
}

class QueryCachedSize{
	int lastQuery = 0;
	int size = 0;

	public QueryCachedSize(int size) {
		this.size = size;
		lastQuery = (int)(Calendar.getInstance().getTimeInMillis()/1000);
	}

	boolean isExpired(int maxsize){
		return size < maxsize || // for small sizes we can always recalculate
				(int)(Calendar.getInstance().getTimeInMillis()/1000) > lastQuery + QSPRConstants.UPDATE_RECORDS_COUNT_TIME;
	}

}	
