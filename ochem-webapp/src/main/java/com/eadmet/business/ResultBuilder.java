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

package com.eadmet.business;

import java.util.List;

import org.hibernate.Criteria;
import org.hibernate.Session;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

public final class ResultBuilder<T>
{

	private final Criteria criteria;
	private final Class<?> rootEntityClass;
	private final Session session; //Hibernate session... "decoupling from Globals.java"
	private PaginationFilter pager;

	private boolean distinct = false;

	public ResultBuilder(Criteria criteria, Class<?> rootEntityClass, Session session) 
	{
		this.criteria = criteria;
		this.rootEntityClass = rootEntityClass;
		this.session = session;
	}

	public ResultBuilder<T> withPager(PaginationFilter pager) 
	{
		this.pager = pager;
		return this;
	}


	public ResultBuilder<T> setDistinct(boolean distinct) 
	{
		this.distinct = distinct;
		return this;
	}


	@SuppressWarnings("unchecked")
	private Criteria getDistinctCriteria()
	{
		criteria.setProjection(Projections.distinct(Projections.id()));
		List<Long> ids = criteria.list();
		Criteria c;
		if (ids.size() > 0)
			c = session.createCriteria(rootEntityClass).add(Restrictions.in("id", ids));
		else
			c = session.createCriteria(rootEntityClass).add(Restrictions.sqlRestriction("1=0"));
		return c;
	}

	@SuppressWarnings("unchecked")
	public List<T> list() 
	{
		Criteria c;
		if (distinct)
			c =  getDistinctCriteria();
		else
			c =  criteria.setProjection(null);

		if (pager != null && pager.pageNum != null && pager.pageSize != null)
		{
			c.setFirstResult((pager.pageNum-1)*pager.pageSize);
			c.setMaxResults(pager.pageSize);
		}

		return c.list();
	}

	public long count() 
	{
		if (distinct)
			return (Long)criteria.setProjection(Projections.countDistinct("id")).list().get(0);
		else
			return (Long)criteria.setProjection(Projections.count("id")).list().get(0);
	}    

	//    public T listOne() 
	//    {
	//    	List<T> list = criteria.list();
	//    	if (list.size() > 0)
	//    		return (T)list.get(0);
	//    	else
	//    		return null;
	//    }
}