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
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.Criteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.dao.Repository;
import qspr.entities.Property;
import qspr.entities.Unit;
import qspr.entities.UnitCategory;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.AccessChecker;
import qspr.util.unitconversion.UnitConversion;
import qspr.workflow.datatypes.DoubleList;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.NumericalValueStandardizer;

@Controller
public class UnitController extends BrowserWrapper
{
	public UnitController()
	{
		sessionRequired = true;
	}

	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		WebList wlist = new WebList();
		Globals.restartAllTransactions(true);

		if (!assertParam("lightweight"))
			Globals.setMarshallingOption(MarshallingOption.UNIT_RECORDS_COUNT);

		Criteria unitCriteria = Globals.session().createCriteria(Unit.class);
		String q = request.getParameter("query");
		String c = request.getParameter("categoryName");
		int psize = 15;
		if(assertParam("query"))
		{
			Disjunction disjunction = Restrictions.disjunction();
			disjunction.add(Restrictions.like("name", "%"+q+"%"));
			disjunction.add(Restrictions.like("description", "%"+q+"%"));
			unitCriteria.add(disjunction);
		}
		if(assertParam("categoryName"))
		{
			unitCriteria.createAlias("category", "category");
			unitCriteria.add(Restrictions.like("category.name", "%"+c+"%"));
		}

		if (assertParam("category"))
		{
			if(!"0".equals(request.getParameter("category")))
			{
				unitCriteria.add(Restrictions.eq("category.id", getLongParam("category")));
				psize = unitCriteria.list().size();
			}			
		} else
			if (assertParam("property"))
			{
				Property property = (Property)Globals.session().get(Property.class, getLongParam("property"));
				if (property != null)
				{
					unitCriteria.add(Restrictions.eq("category.id", property.unitCategory.id));
					psize = unitCriteria.list().size();				
				}
			} else
				psize = 0;

		if (assertParam("pname"))
		{
			unitCriteria
			.createCriteria("category")
			.createCriteria("properties")
			.add(Restrictions.like("shortName", Property.shortName(getParam("pname"))));
			psize = unitCriteria.list().size(); //???
		}

		if (assertParam("id"))
			unitCriteria.add(Restrictions.eq("id", getLongParam("id")));

		if (!assertParam("query"))
			unitCriteria.addOrder(Order.asc("category.id"));
		unitCriteria.addOrder(Order.asc("name"));

		wlist.loadFromCriteria(unitCriteria, getPageNum(), getPageSize(psize));

		//List list = wlist.list;
		Iterator<Object> units = wlist.list.iterator();
		while (units.hasNext())
		{
			Unit unit = (Unit) units.next();
			unit.category.units = null;
			Globals.session().evict(unit);
			//ulist.add(unit);
		}
		Globals.restartAllTransactions(true); // to avoid timeout

		return new WebModel(wlist).getModelAndView();
	}

	public ModelAndView verifyconversion(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Double val0 = NumericalValueStandardizer.getSignificantDigitsDouble(ThreadLocalRandom.current().nextDouble(0.001, 100), NumericalValueStandardizer.SIGNIFICANT_DIGITS + 1);
		Double val1 = UnitConversion.parse(getParam("to-default-conversion"), val0, 100);
		Double val2 = NumericalValueStandardizer.getSignificantDigitsDouble(UnitConversion.parse(getParam("from-default-conversion"), val1, 100),NumericalValueStandardizer.SIGNIFICANT_DIGITS +1);

		List<Double> values = new ArrayList<Double>();
		values.add(val0);
		values.add(val1);
		values.add(val2);

		return new WebModel(new DoubleList(values)).getModelAndView();
	}

	@SuppressWarnings({ "unchecked" })
	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Property p = null;
		if (assertParam("pname"))
		{
			String name = getParam("pname");
			List<Property> l = Globals.session().createCriteria(Property.class).add(Restrictions.eq("shortName", name)).list(); 
			if (l.size() == 0)
				l = Globals.session().createCriteria(Property.class).add(Restrictions.eq("name", name)).list();
			if(l.size() > 0)
				p = l.get(0);

		}
		List <UnitCategory>list = Globals.session().createCriteria(UnitCategory.class).list();
		List<Object> ulist = new ArrayList<Object>();
		for(UnitCategory category: list)
		{
			category.units = null;
			if (p != null && category.equals(p.unitCategory)) { // if we know the property, we should not allow the user to select non-matching units
				category.selected = true;
				ulist.clear();
				Globals.session().evict(category);
				ulist.add(category);
				break;
			}
			Globals.session().evict(category);
			ulist.add(category);
		}
		return new WebModel().setList(ulist).setTemplate("unit-browser").getModelAndView();
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	public ModelAndView edit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Unit unit;

		if (!assertParam("id") || getLongParam("id") < 0)
			unit = new Unit();
		else
			unit = (Unit) Globals.session().get(Unit.class, Long.valueOf(request.getParameter("id")));

		Globals.setMarshallingOption(MarshallingOption.PROPERTY_UNITS);

		List list = Globals.session().createCriteria(UnitCategory.class).list();
		List<Object> ulist = new ArrayList<Object>();
		Iterator<Object> unitcat = list.iterator();
		while(unitcat.hasNext())
		{
			UnitCategory category = (UnitCategory) unitcat.next();
			category.units = null;
			Globals.session().evict(category);
			ulist.add(category);
		}

		return new WebModel(unit).setList(ulist).setRenderMode("popup").setTemplate("unit-edit").getModelAndView();
	}

	public ModelAndView systems(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		return new WebModel().setTemplate("unit-system-browser").getModelAndView();
	}

	public ModelAndView systemslist(HttpServletRequest request, HttpServletResponse response) throws Exception
	{	
		Criteria criteria = Globals.session().createCriteria(UnitCategory.class);
		WebList wlist = new WebList();
		wlist.loadFromCriteria(criteria, getPageNum(), getPageSize(15));
		return new WebModel(wlist).getModelAndView();
	}

	public ModelAndView addNewSystem(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		// FIMXE: Check for superuser privileges
		if (Globals.userSession().user == null)
			throw new UserFriendlyException("Not sufficient privileges to add a new system of units");

		if (!assertParam("defunitname"))
			throw new UserFriendlyException("No default unit specified");

		if (!assertParam("newname"))
			throw new UserFriendlyException("No name specified");

		UnitCategory category = new UnitCategory();
		category.name = getParam("newname");
		Globals.session().save(category);

		Unit unit = new Unit();
		unit.category = category;
		unit.introducer = unit.owner = Globals.userSession().user;
		unit.setName(getParam("defunitname"));
		unit.isdefault = 1;
		Globals.session().save(unit);

		return new WebModel().getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		if(assertParam("action"))
		{
			if (getParam("action").equals("verify"))
				return verifyconversion(request, response);

			boolean moderator = AccessChecker.isModerator(Globals.userSession().user);

			if(!moderator)
				throw new UserFriendlyException("Only moderators can add/edit units. If you need a new one, contact: " + QSPRConstants.INFOEMAIL);

			if (getParam("action").equals("addnewsystem"))
				return addNewSystem(request, response);

			Unit unit;
			if(assertParam("id"))
				unit = (Unit) Globals.session().get(Unit.class, getLongParam("id"));
			else 
				unit = new Unit();

			if (request.getParameter("action").equals("delete"))
				Globals.session().delete(unit);
			else 
			{
				if(unit.id != null)
					unit.doCheckRights();

				unit.owner = Globals.userSession().user;
				if(unit.introducer == null)
					unit.introducer = Globals.userSession().user;
				String name =  request.getParameter("name");
				String description =  request.getParameter("description");
				UnitCategory category = (UnitCategory) Globals.session().get(UnitCategory.class, getLongParam("category"));

				if(!assertParam("id"))
				{
					Unit nunit = Repository.unit.get(name, category);
					if(nunit != null)
						throw new UserFriendlyException("unit already exist in our database. Please perform search before introducing new unit");
				}

				unit.setName(name);
				unit.description =  description;
				unit.category = category;

				if (assertParam("to-default-conversion"))
					unit.toDefaultConversion = getParam("to-default-conversion");
				if (assertParam("from-default-conversion"))
					unit.fromDefaultConversion = getParam("from-default-conversion");
				Globals.session().save(unit);
			}
		}

		return new WebModel().getModelAndView();
	}
}


