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

import java.io.StringReader;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.bind.JAXBException;

import org.hibernate.Criteria;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.entities.Alert;
import qspr.entities.Attachment.AttachmentType;
import qspr.entities.ModelConfigurationTemplate;
import qspr.export.ExportableModel;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.AccessChecker;

@Controller
public class ModelTemplateController extends BrowserWrapper
{
	public ModelTemplateController()
	{
		sessionRequired = true;
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response)
	{
		return new WebModel().setTemplate("model-templates/model-templates-browser").getModelAndView();
	}

	public ModelAndView list(HttpServletRequest request, HttpServletResponse response)
	{
		Criteria c = Globals.session()
				.createCriteria(ModelConfigurationTemplate.class)
				//.add(Restrictions.eq("type", TemplateType.MODEL))
				;
		WebList wl = new WebList();
		wl.loadFromCriteria(c, getPageNum(), getPageSize(15));
		return new BrowserModel()
		.setObject(wl)
		.getModelAndView();
	}

	public ModelAndView edit(HttpServletRequest request, HttpServletResponse response)
	{
		ModelConfigurationTemplate mcTemplate = ModelConfigurationTemplate.getById(getIntParam("id"));
		Globals.setMarshallingOption(MarshallingOption.MODELTEMPLATE_FULLXML);
		return new WebModel(mcTemplate).setTemplate("model-templates/model-templates-edit").getModelAndView();
	}

	public ModelAndView delete(HttpServletRequest request, HttpServletResponse response)
	{
		ModelConfigurationTemplate mcTemplate = ModelConfigurationTemplate.getById(getIntParam("id"));
		mcTemplate.requestModification();
		Globals.session().delete(mcTemplate);

		return new WebModel().getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response)
	{
		if (getParam("action").equals("delete"))
			return delete(request, response);
		return null;
	}

	@SuppressWarnings("unchecked")
	public ModelAndView save(HttpServletRequest request, HttpServletResponse response) throws JAXBException
	{

		ModelConfigurationTemplate mcTemplate = ModelConfigurationTemplate.getById(getIntParam("id"));

		if (!Globals.userSession().user.equals(mcTemplate.introducer))
			AccessChecker.requestModeratorPrivileges();

		ExportableModel model = (ExportableModel) Globals.jaxbContext.createUnmarshaller().unmarshal(new StringReader(getParam("xml")));
		mcTemplate.configuration.setObject(model, AttachmentType.MARSHALABLE);
		mcTemplate.name = getParam("name");
		mcTemplate.updateHash();
		mcTemplate.checkConflicts();

		Globals.session().saveOrUpdate(mcTemplate);

		return redirect("modeltemplate/edit.do?id=" + mcTemplate.id + "&saved=1");
	}

	public ModelAndView getXml(HttpServletRequest request, HttpServletResponse response)
	{
		ModelConfigurationTemplate mcTemplate = ModelConfigurationTemplate.getById(getIntParam("id"));

		return new WebModel(new Alert(mcTemplate.configurationXml)).getModelAndView();
	}

	public ModelAndView setXml(HttpServletRequest request, HttpServletResponse response)
	{
		ModelConfigurationTemplate mcTemplate = ModelConfigurationTemplate.getById(getIntParam("id"));

		//		ExportableModel model = (ExportableModel) mcTemplate.configuration.getObject();
		return new WebModel(new Alert(mcTemplate.configurationXml)).getModelAndView();
	}
}
