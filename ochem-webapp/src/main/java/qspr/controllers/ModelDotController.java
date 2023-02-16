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
import qspr.business.WebFilters;
import qspr.dao.Repository;
import qspr.entities.Alert;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.AccessChecker;

import com.eadmet.business.ModelDotService;
import com.eadmet.business.PaginationFilter;
import com.eadmet.exceptions.UserFriendlyException;

@Controller
public class ModelDotController extends BrowserWrapper
{
	public ModelDotController()
	{
		sessionRequired = true;
	}

	ModelDotService service = new ModelDotService();

	public Model getRequestedModel()
	{
		if (getLongParam("model_id") != null)
			return Repository.model.getById(getLongParam("model_id"));
		else if (assertParam("modelid"))
			return Repository.model.getById(getLongParam("modelid"));
		else
			return (Model) Globals.getSessionAttribute(SessionVariable.MODEL);
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		return new WebModel()
				.setTemplate("model/dot")
				.setRenderMode("popup")
				.getModelAndView();
	}

	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Thread.currentThread().setContextClassLoader(ModelDotController.class.getClassLoader());
		Model model = getRequestedModel();
		long statNum = getLongParam("statnum");
		long epId = getLongParam("ep_id");
		long modelMappingId = getLongParam("mm_id");
		List<ExperimentalProperty> eps = service.getModelDotList(model, statNum, epId, modelMappingId);

		if(!model.approved) {
			ExperimentalProperty ep = Repository.record.getRecord(epId);
			AccessChecker.requestViewingPermission(ep);
		}

		WebList wl = new WebList();
		wl.loadFromList(eps);

		return new WebModel(wl).getModelAndView();		
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public ModelAndView descriptors(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Thread.currentThread().setContextClassLoader(ModelDotController.class.getClassLoader());

		Model model = getRequestedModel();

		List<String> descriptorList = service.getDescriptorList(model, getLongParam("mm_id"), getLongParam("ep_id"));

		WebList wl = new WebList();
		wl.size = descriptorList.size();
		wl.pageNum = getPageNum();
		wl.pageSize = getPageSize(5);
		if ((wl.pageNum-1)*wl.pageSize > descriptorList.size())
			wl.pageNum = 1;
		descriptorList.subList(0, (wl.pageNum-1)*wl.pageSize).clear();
		descriptorList.subList(Math.min(wl.pageSize, descriptorList.size()), descriptorList.size()).clear();
		wl.nativeObjectsList = (List)descriptorList;

		return new WebModel(wl).getModelAndView();		
	}

	public ModelAndView errorshow(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Model model = getRequestedModel();
		return new WebModel(model)
				.setTemplate("model/errors")
				.setRenderMode("popup")
				.addParam("recalculated", "" + assertParam("recalculated"))
				.getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		String action = getParam("action");
		boolean exclude = false;
		if("exclude".equals(action))
			exclude = true;
		Model model = getRequestedModel();
		if (!model.trainingSet.user.equals(Globals.userSession().user))
			throw new UserFriendlyException("You are not authorized to edit basket of user " + model.trainingSet.user.login);

		// We can have multiple IDs separated by comma
		String[] sIds = request.getParameter("id").split(",");
		List<Long> ids = new ArrayList<Long>();
		for (int i = 0; i < sIds.length; i++)
			ids.add(Long.valueOf(sIds[i]));

		long affectedNum = model.trainingSet.excludeOrIncludeEntries(ids, exclude);
		String msg = "" + affectedNum + " records affected";
		return new WebModel(new Alert(msg)).getModelAndView();
	}


	public ModelAndView errorlist(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Thread.currentThread().setContextClassLoader(ModelController.class.getClassLoader());

		PaginationFilter pager = new PaginationFilter(getPageNum(), getPageSize(5));
		Model model = getRequestedModel();

		List<ExperimentalProperty> errors = service.getErrorList(model, pager, assertParam("recalculated"), getParam("errorTitle"));

		WebList webList = new WebList();
		webList.loadFromList(errors);
		webList.pageNum = pager.pageNum;
		webList.pageSize = pager.pageSize;
		webList.size = Long.valueOf(pager.totalSize).intValue();

		WebFilters wf = new WebFilters();
		wf.addFilter("id", model.id.toString(), "");
		if (assertParam("recalculated"))
			wf.addFilter("recalculated", "1", "");
		if (assertParam("errorTitle"))
			wf.addFilter("errorTitle", getParam("errorTitle"), "");

		return new BrowserModel().setFilters(wf).setObject(webList).getModelAndView();
	}

}
