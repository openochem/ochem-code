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

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.Query;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.entities.Journal;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.JournalUtility;

@Controller
public class JournalController extends BrowserWrapper
{
	public JournalController()
	{
		sessionRequired = true;
	}

	public ModelAndView edit(HttpServletRequest request, HttpServletResponse response) 
			throws Exception
	{
		Journal journal;
		if (assertParam("issn"))
			journal = Journal.getByISSN(getParam("issn"));
		else
			journal = (Journal) Globals.session().get(Journal.class, getLongParam("id"));

		if (journal == null)
		{
			journal = new Journal();
			journal.id = Long.valueOf(-1);
			journal.publish_date = null;
		}

		return new WebModel(journal).setTemplate("journal-edit").getModelAndView();
	}


	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) 
			throws Exception
	{
		return new WebModel().setTemplate("journal-browser").getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) 
			throws Exception
	{
		Journal journal = null;
		//get the journal
		if (!request.getParameter("id").equals("-1"))
			journal = (Journal) Globals.session().get(Journal.class, getLongParam("id"));
		else if (assertParam("n-issn"))
			journal = Journal.getByISSN(request.getParameter("n-issn"));
		else if (assertParam("n-title"))
			journal = Journal.getByTitle(request.getParameter("n-title"));

		if (journal == null)
			journal = new Journal();

		if (journal.id != null && ! request.getParameter("action").equals("reload"))
			journal.doCheckRights();

		// delete journal
		if (request.getParameter("action").equals("delete"))
		{
			Globals.session().delete(journal);
		}
		// reload journal
		else if (request.getParameter("action").equals("reload"))
		{
			Globals.session().evict(journal);
			journal.update(JournalUtility.fetchJournal(request.getParameter("n-issn")));

			// Ownership and rights		
			if (Globals.userSession().user != null)
			{
				journal.owner = Globals.userSession().user;
				if (journal.introducer == null)
					journal.introducer = Globals.userSession().user;
			}
		}
		// save journal
		else
		{
			journal.setISSN(request.getParameter("n-issn"));
			journal.setTitle(request.getParameter("n-title"));
			journal.setAbbreviation(request.getParameter("n-abbreviation"));
			journal.publisher = request.getParameter("n-publisher");
			journal.setLink(request.getParameter("n-link"));

			// Ownership and rights		
			if (Globals.userSession().user != null)
			{
				journal.owner = Globals.userSession().user;
				if (journal.introducer == null)
					journal.introducer = Globals.userSession().user;
			}
			Globals.session().saveOrUpdate(journal);
		}

		return new WebModel(journal).getModelAndView();
	}

	public ModelAndView list(HttpServletRequest request, HttpServletResponse response)
	{
		StringBuilder queryText = new StringBuilder("from Journal j where 1 = 1");

		if (assertParam("query"))
			queryText.append(" and (j.title like '%"+getParam("query").toLowerCase()+"%')");

		Query query = Globals.session().createQuery(queryText.toString());

		WebList list = new WebList();
		list.loadFromQuery(query, getPageNum(), getPageSize(10));

		return new BrowserModel().setFilters(formFilters(request)).setObject(list).getModelAndView();
	}
}
