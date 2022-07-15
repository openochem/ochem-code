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

import java.io.IOException;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.entities.Action;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;

@Controller
public class HistoryController extends BrowserWrapper 
{
	public HistoryController()
	{
		sessionRequired = true;
	}
	
	public ModelAndView show(HttpServletRequest req, HttpServletResponse res) throws IOException
	{
		return new WebModel().getModelAndView();
	}
	
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res)
	{
		Criteria criteria = Globals.session().createCriteria(Action.class);
		if (assertParam("record"))
		{
			criteria.add(Restrictions.eq("primaryKey", getLongParam("record")));
			criteria.add(Restrictions.eq("entity", "ExperimentalProperty"));
		}
		criteria.addOrder(Order.desc("id"));
		WebList webList = new WebList();
		webList.loadFromCriteria(criteria, getPageNum(), getPageSize(30));
		
		return new WebModel(webList).getModelAndView();
	}

	public ModelAndView showuser(HttpServletRequest req, HttpServletResponse res) throws IOException
	{
		return new WebModel().setTemplate("userhistory").getModelAndView();
	}
	
	public ModelAndView listuser(HttpServletRequest req, HttpServletResponse res)
	{
		Criteria criteria = Globals.session().createCriteria(Action.class);

		if (assertParam("query"))
			criteria.add(Restrictions.eq("primaryKey", getLongParam("query")));
		
		criteria.createCriteria("useraction")
			.add(Restrictions.eq("user", Globals.userSession().user));
		
		criteria.addOrder(Order.desc("id"));
		WebList webList = new WebList();
		webList.loadFromCriteria(criteria, getPageNum(), getPageSize(30));
		
		return new WebModel(webList).getModelAndView();
	}

}
