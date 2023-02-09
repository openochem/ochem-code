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
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;
import org.hibernate.type.LongType;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.VirtualHostConfiguration;
import qspr.entities.Model;
import qspr.entities.Property;
import qspr.entities.User;
import qspr.frontend.WebModel;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.TimeUtils;

@Controller
public class HomeController extends ControllerWrapper 
{
	/**
	 * Cache for fast loading of the cloud tag
	 */
	private static List<Property> publicProperties = new ArrayList<Property>();
	/**
	 * Cache for counts
	 */
	private static Long totalRecords = null;
	private static Long totalProperties = null;
	private static Long totalArticles = null;
	private static Long countTime = null;

	public ModelAndView index(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException
	{
		if (VirtualHostConfiguration.getHomePage() != null)
			return redirect(VirtualHostConfiguration.getHomePage() + "");
		else
		{
			request.getRequestDispatcher("/index.html").forward(request, response);
			return null;
		}
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) 
			throws Exception 
	{
		if (VirtualHostConfiguration.getHomePage() != null)
			return redirect(VirtualHostConfiguration.getHomePage());

		WebModel wm = new WebModel();

		if (totalRecords == null || (System.nanoTime() - countTime > 6 * 60 * 60 * 1E9)) // Six hours cache invalidation
		{
			countTime = System.nanoTime();
			String where = " from ExperimentalProperty where rights = " + Globals.RIGHTS_FREELY_AVAILABLE + " and property_id !=" + Property.getDummyId();
			totalRecords = Long.valueOf(Globals.session().createSQLQuery("select count(*)" + where).uniqueResult().toString());
			totalProperties = Long.valueOf(Globals.session().createSQLQuery("select count(distinct property_id)" + where).uniqueResult().toString());
			totalArticles = Long.valueOf(Globals.session().createSQLQuery("select count(distinct article_id)" + where).uniqueResult().toString());
		}

		wm.addParam("total-records", "" + totalRecords);
		wm.addParam("total-properties", "" + totalProperties);
		wm.addParam("total-articles", "" + totalArticles);
		wm.addParam("minimum-count", "" + QSPRConstants.MIN_PROPERTY_COUNT);
		return wm.setList(getActiveUsers()).addObjects(getPublishedModels()).setTemplate("home").getModelAndView();
	}

	@SuppressWarnings("unchecked")
	public ModelAndView properties(HttpServletRequest request, HttpServletResponse response)
	{
		synchronized (publicProperties)
		{
			if (publicProperties.isEmpty())
			{
				List<Object[]> rows = Globals.session().createSQLQuery("select property_id, name, count(*) c from ExperimentalProperty join Property using(property_id) " + 
						" where ExperimentalProperty.rights = " + Globals.RIGHTS_FREELY_AVAILABLE + " and property_id != " + Property.getDummyId() + " group by ExperimentalProperty.property_id having c >= " + QSPRConstants.MIN_PROPERTY_COUNT )
						.list();
				for (Object[] row : rows) 
				{
					if (row[0] == null)
						continue;
					Property p = new Property();
					p.id = ((Number) row[0]).longValue();
					p.setName((String) row[1]);
					p.count = ((Number) row[2]).longValue();
					publicProperties.add(p);
				}
			}
		}
		return new WebModel().setList(publicProperties).getModelAndView();
	}

	@SuppressWarnings("unchecked")
	private List<User> getActiveUsers()
	{
		String loginsToHide = "", ipToHide = "";

		for(String login : User.getDevelopers())
			loginsToHide += " and BINARY login != '" + login + "'";

		String ips[] = {"0:0:0:0:0:0:0:1"};

		for(String ip: ips)
			ipToHide += " and ip_address not like '"+ ip+"'";

		String query = "select max(activity_time) m, user_id from Session natural left join User where "
				+ "(user_id is not null " 
				+ loginsToHide
				+ ipToHide
				+ ") group by user_id order by m desc limit 6;\n";

		//		System.out.println(query);

		List<Object[]> rows = Globals.session().createSQLQuery(query).list();

		List<User> users = new ArrayList<User>();
		for (Object[] row : rows) {
			User user = (User) Globals.session().get(User.getCurrentClass(), ((Number)row[1]).longValue());
			user.latestActivity = TimeUtils.ago((Timestamp) row[0]);
			users.add(user);
		}

		return users;
	}



	@SuppressWarnings("unchecked")
	private List<Model> getPublishedModels()
	{
		List<Long> modelIds = Globals.session().createSQLQuery("select max(model_id) mx from Model natural left join Session where published=1 and approved=1 group by user_id").addScalar("mx", LongType.INSTANCE).list();
		List<Model> models = Globals.session().createCriteria(Model.class).add(Restrictions.in("id", modelIds)).addOrder(Order.desc("lastModification")).list();

		Set<Model> depdupeModelIds = new LinkedHashSet<Model>(models);
		models.clear();
		models.addAll(depdupeModelIds);

		return models;
	}
}
