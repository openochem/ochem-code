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
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;

import org.hibernate.Query;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.entities.Author;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;

@Controller
public class AuthorController extends BrowserWrapper 
{
	Author author;

	public AuthorController()
	{
		sessionRequired = true;
	}
	
	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) 
	throws Exception
	{
		return new WebModel().setTemplate("author-browser").getModelAndView();
	}
	
	public ModelAndView edit(HttpServletRequest request, HttpServletResponse response) 
	throws Exception
	{
		Author auth = (Author) Globals.session().get(Author.class, getLongParam("id"));
		if (auth == null)
		{
			auth = new Author();
			auth.id = Long.valueOf(-1);
		}
		return new WebModel(auth).setTemplate("author-edit").getModelAndView();
	}
	
	public ModelAndView action(HttpServletRequest req, HttpServletResponse res)
	{
		Author auth = (Author) Globals.session().get(Author.class, getLongParam("id"));
		if (auth == null)
		{
			auth = new Author();
			auth.id = Long.valueOf(-1);
		}
		if (req.getParameter("action").equals("delete"))
		{
			Globals.session().delete(auth);	
			return new WebModel().getModelAndView();
		} else
		{
			if (assertParam("n-firstname"))
				auth.firstName = req.getParameter("n-firstname");
			
			if (assertParam("n-lastname"))
				auth.lastName = req.getParameter("n-lastname");
			
			if (assertParam("n-initials"))
				auth.initials = req.getParameter("n-initials");
	
			Globals.session().save(auth);
			return new WebModel(auth).getModelAndView();
		}
		
	}
	
	public ModelAndView list(HttpServletRequest request, HttpServletResponse response)
	{
		String q = request.getParameter("query");
		if (q == null)
			q = "";
		
		Query query = Globals.session().createQuery("from Author p where lower(p.lastName) like '%"+q+"%' ");
		WebList list = new WebList();
		list.loadFromQuery(query, getPageNum(), getPageSize(15));
		return new WebModel(list).getModelAndView();
	}
}

@XmlType(name = "author-list")
class AuthorList
{
	@XmlElement(name = "author")
	List<Author> authors = new ArrayList<Author>();
}
