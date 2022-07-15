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

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.entities.CalculatedDescriptor;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.modelling.configurators.BasicModelConfigurator;
import qspr.modelling.configurators.DescriptorsConfigurator;

@Controller
public class DescriptorsController extends BrowserWrapper 
{
	public DescriptorsController()
	{
		sessionRequired = true;
	}
	
	private List<CalculatedDescriptor> getWorkingSet()
	{
		List<CalculatedDescriptor> filtered = new ArrayList<CalculatedDescriptor>();
		
		BasicModelConfigurator object = (BasicModelConfigurator) Globals.getSessionAttribute(SessionVariable.MODEL_CONFIGURATOR);
		if (object == null || !(object instanceof DescriptorsConfigurator))
			return filtered;
		
		List<CalculatedDescriptor> original = ((DescriptorsConfigurator) Globals.getSessionAttribute(SessionVariable.MODEL_CONFIGURATOR)).manualDescriptorList;
		
		for (CalculatedDescriptor d : original)
		{
			if (assertParam("selected") && !d.selected)
				continue;
			if (assertParam("name") && !d.name.matches(".*"+getParam("name")+".*"))
				continue;
			if (assertParam("id") && !d.id.equals(Long.valueOf(getParam("id"))))
				continue;			
			filtered.add(d);
		}
		return filtered;
				
	}
	
	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		return new WebModel().setTemplate("modeller/configurators/select-descriptors").getModelAndView();
	}
	
	public ModelAndView action(HttpServletRequest request, HttpServletResponse response)
	{
		String action = getParam("action");
		List<CalculatedDescriptor> filtered = getWorkingSet();
		for (CalculatedDescriptor d : filtered) 
		{
			if (action.equals("select"))
				d.selected = true;
			else if (action.equals("unselect"))
				d.selected = false;
			else if (action.equals("toggle"))
				d.selected = !d.selected; 
				
		}
		return new WebModel().getModelAndView();
	}
	
	public ModelAndView list(HttpServletRequest request, HttpServletResponse response)
	{
		List<CalculatedDescriptor> filtered = getWorkingSet();
		WebList wl = new WebList();
		wl.loadFromList(filtered, getPageNum(), getPageSize(15));
		return new WebModel(wl).getModelAndView();
		
//		Criteria criteria = Globals.session().createCriteria(CalculatedDescriptor.class);
//		criteria.add(Restrictions.eq("session", Globals.userSession()));
//		
//		if (assertParam("selected"))
//		{
//			criteria.add(Restrictions.or(Restrictions.eq("selected", true), Restrictions.eq("status", 1)));
//		}
//		
//		if (assertParam("name"))
//			criteria.add(Restrictions.like("name", "%"+getParam("name")+"%"));
//		
//		WebList wl = new WebList();
//		wl.loadFromCriteria(criteria, getPageNum(), getPageSize(15));
//		
//		return new WebModel(wl).getModelAndView();
	}
}

