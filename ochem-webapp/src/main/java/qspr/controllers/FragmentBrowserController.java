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

import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.business.WebFilters;
import qspr.entities.Molecule;
import qspr.fragmententities.Mapping1Fragment;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;

@Controller
public class FragmentBrowserController extends BrowserWrapper 
{
	private static transient final Logger logger = LogManager.getLogger(FragmentBrowserController.class);

	public FragmentBrowserController()
	{
		sessionRequired = true;
	}

	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		WebFilters filters = formFilters(request);
		WebList list = new WebList();
		list.column = "mapping1_id";
		Criteria mainCriteria = createMainCriteria(filters, Globals.alternateSession().createCriteria(Mapping1Fragment.class, "mf"));
		list.useEntity(Mapping1Fragment.class).loadFromCriteria(mainCriteria, getPageNum(), getPageSize(5));
		return new BrowserModel().setFilters(filters).setObject(list).getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception{

		return new WebModel().getModelAndView();
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Long id = getLongParam("id");
		Criteria c = Globals.session().createCriteria(Molecule.class)
				.createAlias("mapping1","mp1")
				.add(Restrictions.eq("mp1.id",id))
				.setMaxResults(1);
		@SuppressWarnings("unchecked")
		List<Molecule> mlist = c.list();
		Molecule m = null;
		if (mlist.size() > 0)
			m = mlist.get(0);

		return new WebModel().setObject(m).setTemplate("fragment-browser").getModelAndView();
	}



	private Criteria createMainCriteria(WebFilters filters, Criteria criteria )
	{
		if (filters.has("id"))
		{
			// ID has been explicitly specified
			// No need to look for other filters

			String[] ids = ThreadScope.get().localRequest.getParameterValues("id");
			if (ids == null)
				ids = new String[]{filters.get("id")};
			if (ids.length == 1)
			{	logger.info(filters.get("id"));
			criteria.add(Restrictions.eq("mapping1_id", filters.getLong("id")));
			//criteria.setMaxResults(getPageSize(5)).setFirstResult((getPageNum() -1)*getPageSize(5)).list();

			//				logger.info("List query executed: " + stats.current());
			criteria.createAlias("fragment", "frag");
			//				criteria.add(Restrictions.eq("mf.fragment_id", "frag.id"));
			}

			Order order = (filters.has("order")) ? Order.asc("frag.size") : Order.desc("frag.size");	
			if (filters.has("sortby") &&filters.get("sortby").equals("number"))
			{
				order = (filters.has("order")) ? Order.asc("mf.fragment_count") : Order.desc("mf.fragment_count");
			}
			criteria.addOrder(order);

		}
		return criteria;

	}
}
