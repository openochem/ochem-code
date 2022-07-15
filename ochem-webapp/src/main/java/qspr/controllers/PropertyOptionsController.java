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

import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.mailer.Mailer;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.AccessChecker;

/**
 * Controller to manage property options
 * @author midnighter
 */
@Controller
public class PropertyOptionsController extends BrowserWrapper
{
	public PropertyOptionsController()
	{
		sessionRequired = true;
	}

	@Override
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		Criteria c = Globals.session().createCriteria(PropertyOption.class).add(Restrictions.eq("property", Property.getById(getLongParam("property"))));
		if (assertParam("name"))
			c.add(Restrictions.like("name", "%" + getParam("name") + "%"));

		c.addOrder(Order.asc("name"));
		WebList wl = new WebList();
		wl.loadFromCriteria(c, getPageNum(), getPageSize(50));
		return new BrowserModel().setObject(wl).getModelAndView();
	}

	public ModelAndView show(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		Property p = Property.getById(getLongParam("property"));
		return new WebModel(p).setTemplate("browsers/options-browser").getModelAndView();
	}

	public ModelAndView add(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		Property p = Property.getById(getLongParam("property"));
		AccessChecker.requestModificationPermission(p);

		// Can be several options splitted with a newline sign
		String[] names = getParam("name").split("[\n\r]+");
		for (String name : names)
		{
			name = name.trim();
			if (!"".equals(name))
			{
				if (p.getOption(name) == null)
				{
					PropertyOption po = Repository.option.getPropertyOptionByName(name, p.id, true, false);
					Globals.session().save(po);
				}
			}
		}

		return new WebModel().getModelAndView();
	}

	public ModelAndView edit(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		long optionId = getLongParam("id");
		PropertyOption po = Repository.option.getPropertyOptionById(optionId);
		AccessChecker.requestModificationPermission(po.property);
		String name = getParam("name").trim();
		System.out.println(optionId + " " +name);
		if(!po.name.equals(name)){ // only a new name
			PropertyOption pold =  Repository.option.getPropertyOptionByName(name, po.property.id, false, false);
			if(pold != null)
				throw new UserFriendlyException("Option with such name alraedy exist");
			Mailer.notifyAdmins("Option renamed for property " + po.property.getName(), "Old : " + po.name + " to new: " + name + " by " + Globals.userSession().user.login);
			po.name = name;
			Globals.session().saveOrUpdate(po);
		}
		return new WebModel().getModelAndView();
	}

	public ModelAndView delete(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		PropertyOption po = (PropertyOption) Globals.session().get(PropertyOption.class, getLongParam("id"));
		AccessChecker.requestModificationPermission(po.property);
		Globals.session().delete(po);

		return new WebModel().getModelAndView();
	}

}
