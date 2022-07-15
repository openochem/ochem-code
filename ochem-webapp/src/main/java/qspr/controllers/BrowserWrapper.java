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

import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;
import org.springframework.web.servlet.ModelAndView;

import qspr.ThreadScope;
import qspr.business.WebFilters;

import com.eadmet.business.PaginationFilter;

public abstract class BrowserWrapper extends ControllerWrapper 
{
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res) throws Exception {
		return null;
	}

	public int getPageNum()
	{
		if (getParam("pagenum") != null)
			return getIntParam("pagenum");
		else
			return 1;
	}

	public int getPageSize(int defaultValue)
	{
		if (assertParam("pagesize"))
			return Math.min(getIntParam("pagesize"), 500);
		else
			return defaultValue;
	}

	@SuppressWarnings("unchecked")
	protected WebFilters formFilters(HttpServletRequest request)
	{
		WebFilters filters = new WebFilters();
		Enumeration<String> params = request.getParameterNames();
		while (params.hasMoreElements()) 
		{
			String key = params.nextElement();
			String value = request.getParameter(key).trim();
			if ("out".equals(key) || "render-mode".equals(key))
				continue;
			formFilter(key, value, filters);
		}		
		return filters;
	} 

	protected void formFilter(String key, String value, WebFilters filters)
	{
		filters.addFilter(key, value, null);	
	}

	protected void filterById(Criteria criteria)
	{
		String[] ids = ThreadScope.get().localRequest.getParameterValues("id");
		List<Long> idsList = new ArrayList<Long>();

		for (String id : ids)
			if (id.contains(","))
			{
				// A comma-separated list of IDs
				String[] commaIds = ids[0].split(",");

				for (String idInner : commaIds) 
					idsList.add(Long.valueOf(idInner));
			}
			else
				idsList.add(Long.valueOf(id));
		if (idsList.size() > 1)
			criteria.add(Restrictions.in("id", idsList));
		else
			criteria.add(Restrictions.eq("id", idsList.get(0)));
	}

	public static PaginationFilter getPaginationFilter(HttpServletRequest r)
	{
		PaginationFilter pf = new PaginationFilter();
		if (r.getParameter("pageNum") != null)
			pf.pageNum = Integer.valueOf((String)r.getParameter("pageNum"));
		else if (r.getParameter("pagenum") != null)
			pf.pageNum = Integer.valueOf((String)r.getParameter("pagenum"));
		else 
			pf.pageNum = 1;

		if (r.getParameter("pageSize") != null)
			pf.pageSize = Integer.valueOf((String)r.getParameter("pageSize"));
		else if (r.getParameter("pagesize") != null)
			pf.pageSize = Integer.valueOf((String)r.getParameter("pagesize"));
		else
			pf.pageSize = 15;

		return pf;
	}

}
