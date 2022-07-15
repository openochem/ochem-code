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

import java.io.File;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.OCHEMConfiguration;
import qspr.entities.ExportAction;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.metaserver.protocol.NoSQLReference;
import qspr.metaserver.transport.NoSqlTransport;
import qspr.util.ExportThread;
import qspr.util.WrapperThread;

import com.eadmet.utils.FileUtils;

@Controller
public class ExportController extends BrowserWrapper 
{

	//	private static transient final Logger logger = Logger.getLogger(ExportController.class);

	public ExportController()
	{
		sessionRequired = true;
	}

	public ExportThread findThread(ExportAction ea)
	{
		for (WrapperThread wt : WrapperThread.runningThreads)
			if (wt instanceof ExportThread)
			{
				ExportAction eaction = ((ExportThread)wt).eAction;
				if (eaction != null && eaction.id.equals(ea.id))
					return (ExportThread)wt;
			}
		return null;
	}

	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		WebList wlist = new WebList();

		Criteria exportHistoryCriteria = Globals.session().createCriteria(ExportAction.class).createAlias("session", "s");
		exportHistoryCriteria.add(Restrictions.eq("s.user", Globals.userSession().user)).add(Restrictions.isNotNull("fileName")).add(Restrictions.isNotNull("status")).addOrder(Order.desc("id"));
		wlist.loadFromCriteria(exportHistoryCriteria, getPageNum(), getPageSize(50));
		for (Object obj : wlist.list)
		{
			ExportAction ea = (ExportAction)obj;
			if (!ea.failed && (!"Finished".equals(ea.status) && ea.status != null))
			{
				ExportThread et = findThread(ea);
				if (et != null)
				{
					ea.running = true;
					ea.status = et.getStatus();
				} else
				{
					ea.failed = true;
					ea.running = false;
					ea.status = "Export crashed with unrecoverable error";
				}
			}
		}
		return new WebModel(wlist).getModelAndView();
	}

	public ModelAndView confirm(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ExportAction ea = (ExportAction)Globals.session().get(ExportAction.class, getLongParam("id"));
		Globals.setMarshallingOption(MarshallingOption.EXPORTACTION_BONUSES);
		//Transaction, free bonus points used, etc

		ea.fileDownloaded = true;

		return WebModel.redirect("export/download.do?id="+ea.id).getModelAndView();
	}

	public ModelAndView download(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ExportAction ea = (ExportAction)Globals.session().get(ExportAction.class, getLongParam("id"));
		if (OCHEMConfiguration.inhouseInstallation)
			ea.fileDownloaded = true;
		if (ea.fileDownloaded) //He already accepted the transaction
		{
			String s = ea.getFullFilePath();
			if (!new File(s).exists())
			{
				//Get from MongoDB and put to the folder				
				NoSQLReference ref = new NoSQLReference(null, "exports", "exports");
				byte[] file = NoSqlTransport.getDataSafely(ref);
				FileUtils.saveBytesToFile(file, s);
			}
		}
		Globals.setMarshallingOption(MarshallingOption.EXPORTACTION_BONUSES);
		return new WebModel(ea).setTemplate("file-download").getModelAndView();
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		return new WebModel().setTemplate("exporthistory-browser").getModelAndView();
	}
}
