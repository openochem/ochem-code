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

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.export.ExportableColumn;
import qspr.export.ExportableSet;
import qspr.export.ExportableSetConfiguration;
import qspr.frontend.WebModel;
import qspr.util.ExportThread;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.business.DescriptorsStorageService;
import com.eadmet.business.DescriptorsStorageService.CacheUploadResult;
import com.eadmet.exceptions.UserFriendlyException;

@Controller
public class DescriptorsStorageController extends ControllerWrapper
{

	DescriptorsStorageService service = new DescriptorsStorageService();

	public DescriptorsStorageController()
	{
		sessionRequired = true;
	}

	public ModelAndView show(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		if (!Globals.isValidatedUser())
			throw new UserFriendlyException(" Only validated experienced users can use this functionality. Most of users should NEVER upload any descriptors: descriptors i OCHEM are calculated on-line from chemical structures." +
					" Uploaded descriptors canonly be used in very specific situations, e.g. for comparison with other descriptor packages which are not available in the OCHEM. Otherwise is  on-sence to upload them. " + 
					" Absent packages can be added on demand. If you need it, contact as at: " +
					QSPRConstants.INFOEMAIL
					);
		WebModel wm =  new WebModel();
		if (Globals.isSuperUser())
			wm.setList(service.getPublicConfigs());
		else
			wm.setList(service.getPrivateConfigs());
		return wm.setTemplate("descriptorsstorage/overview").getModelAndView();
	}

	public ModelAndView delete(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		service.deleteCache(getParam("config-id"), getParam("user"));
		return new WebModel().getModelAndView();
	}

	public ModelAndView uploadSubmit(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		if (!Globals.isValidatedUser())
			throw new UserFriendlyException(" Only validated experienced users can use this functionality. Most of users should NEVER upload any descriptors: descriptors i OCHEM are calculated on-line from chemical structures." +
					" Uploaded descriptors canonly be used in very specific situations, e.g. for comparison with other descriptor packages which are not available in the OCHEM. Otherwise is  on-sence to upload them. " + 
					" Absent packages can be added on demand. If you need it, contact as at: " +
					QSPRConstants.INFOEMAIL
					);

		CacheUploadResult result = service.descriptorsCacheUpload(Globals.getUploadedFile(), getParam("desc-type"), getParam("desc-conf-xml"));
		// For UI
		WebModel wm = new WebModel();
		if (result.externalIDsUpload)
			wm.addParam("external-mols", "true");

		return wm.setTemplate("descriptorsstorage/upload-success").getModelAndView();
	}

	public ModelAndView export(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		ExportableSet eData = new ExportableSet();
		eData.clearColumns();
		eData.addColumn(ExportableColumn.SMILES);
		eData.addColumn(ExportableColumn.MOLECULEID);
		eData.addColumn(ExportableColumn.EXTERNAL_ID);
		eData.addColumn(ExportableColumn.DESCRIPTORS);

		if (assertParam("submit"))
		{
			String user = assertParam("public") ? null : Globals.userSession().user.login;
			ExportThread eThread = service.getExportThread(getParam("config-id"), user, ExportableSetConfiguration.configureFromDialog(req), getParam("format"));
			eThread.start();
			return redirect("longoperations/operationWaitingScreen.do?operation-id=" + eThread.operationID);
		}
		else
			return new WebModel(eData).setTemplate("export").getModelAndView();
	}
}
