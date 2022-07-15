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

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.entities.Article;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.User;
import qspr.frontend.WebModel;
import qspr.util.AccessChecker;

import com.eadmet.datacube.DataTree;
import com.eadmet.datacube.DataTree.Callback;
import com.eadmet.datacube.Subkonto;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.moderation.AwaitingDataReport;
import com.eadmet.moderation.DataMetrics;
import com.eadmet.moderation.DataReportRequest;
import com.eadmet.moderation.ModeratorService;

@Controller
public class EditorController extends ControllerWrapper 
{
	
	
	public ModelAndView getReportData(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Globals.setMarshallingOption(MarshallingOption.PROPERTY_MODERATOR);
		AwaitingDataReport report;
		
		DataReportRequest reportRequest = new DataReportRequest();
		
		reportRequest.groupingOrder = new ArrayList<Subkonto>();
		String[] groupings = request.getParameter("groupings").split(",");
		for (String grouping : groupings)
		{
			if (grouping.equals("cond-"))
				continue;
			reportRequest.groupingOrder.add(getSubkontoByName(grouping));
			if (grouping.startsWith("cond-"))
				reportRequest.qualitativeCondition = Property.getById(Long.valueOf(grouping.substring(5)));
		}
		
		reportRequest.countPrivateRecords = false;
		if (assertParam("property")) 
		{
			reportRequest.properties.add(Property.getById(getLongParam("property")));
			reportRequest.countPrivateRecords = Globals.isSuperUser();
		}
		if (assertParam("my"))
		{	
			reportRequest.user = Globals.userSession().user;
			reportRequest.countPrivateRecords = true;
		}
		else if (assertParam("introducer"))
		{
			reportRequest.user = User.getByString(getParam("introducer"));
			reportRequest.countPrivateRecords = Globals.isSuperUser();
		}
		else if (!Globals.isSuperUser() && assertParam("moderator"))
			reportRequest.moderator = Globals.userSession().user;
		
		reportRequest.countApprovedRecords = (reportRequest.user != null || reportRequest.moderator != null || !reportRequest.properties.isEmpty());
		
		checkPrivileges(reportRequest);
		
		if (assertParam("publicOnly"))
			reportRequest.countPrivateRecords = false;
		
		report = moderatorService.getAwaitingDataReport(reportRequest);
		
		// Reduce the size of the tree
		report.dataTree.walk(new Callback<DataTree<DataMetrics>>()
		{
			@Override
			public void apply(DataTree<DataMetrics> node)
			{
				if (node.currentGrouping instanceof Article) {
					Article a = (Article) node.currentGrouping;
					Article reduced = new Article();
					reduced.setTitle(a.getTitle());
					reduced.id = a.id;
					node.currentGrouping = reduced;
				}
			}
		});
		
		return new WebModel(report).addObject(reportRequest).getModelAndView();
	}
	
	private Subkonto getSubkontoByName(String name) {
		name = name.toLowerCase();
		if (name.isEmpty())
			throw new UserFriendlyException("Please, select at least one grouping");
		if ("introducer".equals(name))
			return new Subkonto(User.getCurrentClass(), "Introducer");
		else if ("property".equals(name))
			return new Subkonto(Property.class, "Property");
		else if ("article".equals(name))
			return new Subkonto(Article.class, "Article");
		else if ("isprimary".equals(name))
			return new Subkonto(Boolean.class, "Primary", "isPrimary");
		else if (name.startsWith("cond-"))
			return new Subkonto(PropertyOption.class, Property.getById(Long.valueOf(name.substring(5))).getName(), "pv.option");
		else
			throw new UserFriendlyException("Unsupported grouping: " + name);
	}
	
	public EditorController()
	{
		sessionRequired = true;
	}
	
	private void checkPrivileges(DataReportRequest request)
	{
		request.countPrivateRecords = request.user != null && (request.user.equals(Globals.myself()) || Globals.isSuperUser());
		if (Globals.isSuperUser())
			return;
		if (request.user != null && request.user.equals(Globals.userSession().user))
			return;
		if (AccessChecker.isModerator(Globals.userSession().user) && request.moderator != null && request.moderator.equals(Globals.userSession().user))
			return;
		if (request.user != null)
			return;
		if (!request.properties.isEmpty())
			return;
		
		throw new UserFriendlyException("You cannot access this report");
	}
	
	@Autowired
	private ModeratorService moderatorService;
}
