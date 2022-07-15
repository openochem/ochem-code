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
import java.io.OutputStream;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.business.BatchUpoadBrowserFilter;
import qspr.dao.Repository;
import qspr.entities.Alert;
import qspr.entities.User;
import qspr.frontend.MarshalableList;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.batchupload.entityschema.EntitiesRemapping;
import com.eadmet.batchupload.main.BatchUploadProcessor;
import com.eadmet.batchupload.main.MultithreadBatchUploadProcessor;
import com.eadmet.batchupload.main.RecordPreview;
import com.eadmet.batchupload.main.RecordPreview.PreviewRecordStatus;
import com.eadmet.batchupload.main.RecordPreview.PreviewUploadAction;
import com.eadmet.batchupload.main.TransactionBatchUploadEventHandler;
import com.eadmet.batchupload.main.UploadPreview;
import com.eadmet.batchupload.main.UploadedColumnSchema;
import com.eadmet.batchupload.main.UploadedColumnSchema.UploadedColumnType;
import com.eadmet.batchupload.main.UploadedFileSchema;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;

@Controller
public class BatchUpload30Controller extends ControllerWrapper 
{
	public static final String BATCH_PROCESSOR = "batch-processor";
	public static final String BATCH_ERROR = "batch-error";
	public static final String BATCH_PAGE = "batch-page";

	public BatchUpload30Controller()
	{
		sessionRequired = true;
	}

	protected MultithreadBatchUploadProcessor getProcessor()
	{
		MultithreadBatchUploadProcessor processor = (MultithreadBatchUploadProcessor)ThreadScope.get().localRequest.getSession().getAttribute(BATCH_PROCESSOR);
		if (processor == null)
		{
			processor = new MultithreadBatchUploadProcessor();
			processor.eh = new TransactionBatchUploadEventHandler();
			ThreadScope.get().localRequest.getSession().setAttribute(BATCH_PROCESSOR, processor);
		}
		return processor;
	}

	protected void setCurrentPage(String page)
	{
		ThreadScope.get().localRequest.getSession().setAttribute(BATCH_PAGE, page);
	}

	protected String getSensiblePageDefaults(String page)
	{
		BatchUploadProcessor processor = getProcessor();
		if (processor.parser == null)
			return "init";
		return page;
	}

	protected String getCurrentPage()
	{
		String page = (String)ThreadScope.get().localRequest.getSession().getAttribute(BATCH_PAGE);
		if (page == null)
			return "init";
		else
			return getSensiblePageDefaults(page);
	}

	public ModelAndView show(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		if (req.getParameter("nodb") == null && Globals.userSession().user == null)
			return redirect("user/pleaseregister.do");

		return new ModelAndView("redirect:"+getCurrentPage()+".do?render-mode=popup");
	}

	public ModelAndView cancel(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		setCurrentPage(null);
		ThreadScope.get().localRequest.getSession().setAttribute(BATCH_PROCESSOR, null);
		return new ModelAndView("redirect:show.do?render-mode=popup&nodb=true");
	}

	public ModelAndView init(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		setCurrentPage("init");
		BatchUploadProcessor processor = getProcessor();
		try
		{
			String error = (String)ThreadScope.get().localRequest.getSession().getAttribute(BATCH_ERROR);
			return new WebModel(processor.context).addParam("error", error).setTemplate("batch30/init").getModelAndView();
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}

	public ModelAndView init_submit(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		BatchUploadProcessor processor = getProcessor();
		try
		{
			processor.context.allowArticlePubmedSearch = (req.getParameter("allow_pubmed") != null);
			processor.context.allowMoleculePubchemSearch = (req.getParameter("allow_pubchem") != null);
			processor.context.hiddenByDefault = (req.getParameter("hidden") != null);
			File f = Globals.getUploadedFile();
			processor.setFile(f);
			setCurrentPage("remap_columns");
		} catch (Exception e)
		{
			e.printStackTrace();
			setCurrentPage("init");
			ThreadScope.get().localRequest.getSession().setAttribute(BATCH_ERROR, OCHEMUtils.exceptionToString(e));
		}
		return new ModelAndView("redirect:show.do?render-mode=popup&nodb=true");
	}

	public ModelAndView remap_columns(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		setCurrentPage("remap_columns");
		BatchUploadProcessor processor = getProcessor();
		try
		{
			return new WebModel(processor.getFileSchema())
					.addObject(new MarshalableList(UploadedColumnType.getKnownColumns()).setName("known"))
					.setTemplate("batch30/remap_columns").getModelAndView();
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}

	public ModelAndView validate_fileschema(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		BatchUploadProcessor processor = getProcessor();
		try 
		{
			UploadedFileSchema ufs = processor.getFileSchema();
			ufs.setSelectedSheet(getIntParam("sheet"));
			String[] names = req.getParameterValues("name");
			for (int i=0; i<names.length; i++)
			{
				UploadedColumnSchema ucs = ufs.getSelectedSheetSchema().getColumns().get(i);
				ucs.ignore = (ThreadScope.get().localRequest.getParameter("select"+i) == null);
				ucs.name = names[i];
			}

			if (assertParam("adddummy"))
				ufs.getSelectedSheetSchema().addDummyProperty = true;
			else
				ufs.getSelectedSheetSchema().addDummyProperty = false;

			ufs = processor.getFileSchema();
			return new WebModel(ufs).getModelAndView();
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}

	public ModelAndView remap_columns_submit(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		BatchUploadProcessor processor = getProcessor();
		try
		{
			UploadedFileSchema ufs = processor.getFileSchema();

			ufs.setSelectedSheet(getIntParam("sheet"));

			String[] names = req.getParameterValues("name");

			for (int i=0; i<names.length; i++)
			{
				UploadedColumnSchema ucs = ufs.getSelectedSheetSchema().getColumns().get(i);
				ucs.ignore = (ThreadScope.get().localRequest.getParameter("select"+i) == null);
				ucs.name = names[i];
			}

			if (assertParam("adddummy"))
				ufs.getSelectedSheetSchema().addDummyProperty = true;
			else
				ufs.getSelectedSheetSchema().addDummyProperty = false;

			processor.getRemappingSchema();
			setCurrentPage("remap_entities");
			return new ModelAndView("redirect:show.do?render-mode=popup&nodb=true");
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}


	public ModelAndView remap_entities(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		BatchUploadProcessor processor = getProcessor();
		try 
		{
			if (processor.getRemappingSchema() == null)
			{
				return new ModelAndView("redirect:status.do?render-mode=popup&nodb=true");
			} else
			{
				EntitiesRemapping ers = processor.getRemappingSchema();
				return new WebModel(ers).setTemplate("batch30/remap_entities").getModelAndView();
			}
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}

	@SuppressWarnings("unchecked")
	public ModelAndView remap_entities_submit(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		MultithreadBatchUploadProcessor processor = getProcessor();
		try
		{
			EntitiesRemapping ers = processor.getRemappingSchema();
			ers.remapFromParameterMap(req.getParameterMap());
			processor.getRemappingSchema();
			processor.waitToFinish();
			boolean warningIsInvalid = !assertParam("ignorewarnings");
			if (processor.getRemappingSchema().valid(warningIsInvalid))
			{
				processor.getUploadPreview();
				setCurrentPage("browser");
			}
			return new ModelAndView("redirect:show.do?render-mode=popup&nodb=true");
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}

	public ModelAndView browser(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		BatchUploadProcessor processor = getProcessor();
		try 
		{
			if (processor.getUploadPreview() == null)
			{
				return new ModelAndView("redirect:status.do?render-mode=popup&nodb=true");
			} else
			{
				UploadPreview preview = processor.getUploadPreview();
				return new WebModel(preview.getSummary()).setTemplate("batch30/browser").getModelAndView();
			}
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}

	public ModelAndView browser_submit(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		BatchUploadProcessor processor = getProcessor();

		User user = Globals.userSession().user;
		int records = Repository.record.getUnpublishedRecordsCountForUser(user.id);
		int num = 3;

		if(!user.isOCHEMDeveloper() && !user.isSuperUser() && !OCHEMConfiguration.inhouseInstallation) {
			if(user.isValidated() && records > num*QSPRConstants.VALIDATEDUSER_RECORDS_LIMIT)
				throw new UserFriendlyException(" As a validated user you can upload and store maximum " + num*QSPRConstants.VALIDATEDUSER_RECORDS_LIMIT + " non-published records, but you have n=" + records + " hidden records. Review your old/unused data and delete non-required ones or make them publicly available.");
			if(!user.isValidated() && records > num*QSPRConstants.USER_RECORDS_LIMIT)
				throw new UserFriendlyException(" As a non-validated user you can upload and store  maximum " + num*QSPRConstants.USER_RECORDS_LIMIT + " non-published records, but you have n=" + records + " hidden records. Review your old/unused data and delete non-required ones, make them publicly available or ask to validate your account at " +QSPRConstants.INFOEMAIL);
		}

		try 
		{
			processor.upload();
			setCurrentPage("result");
			return new ModelAndView("redirect:show.do?render-mode=popup&nodb=true");
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}			
	}

	public ModelAndView result(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		BatchUploadProcessor processor = getProcessor();
		try
		{
			if (processor.upload() == null)
			{
				return new ModelAndView("redirect:status.do?render-mode=popup&nodb=true");
			} else
			{
				UploadPreview preview = processor.upload();
				return new WebModel(preview.getSummary()).setTemplate("batch30/result").getModelAndView();
			}
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}

	public ModelAndView report(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		BatchUploadProcessor processor = getProcessor();
		try
		{
			byte[] data = processor.getUploadReport();
			res.setContentType("application/xls");
			res.setHeader("Content-Disposition", "attachment; filename="+processor.parser.sparser.fileName+"-upload-report.xls");
			OutputStream os = res.getOutputStream();
			os.write(data);
			os.flush();
			os.close();
			return null;
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public ModelAndView browser_list(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		BatchUploadProcessor processor = getProcessor();
		try
		{
			Integer pageNum = (getParam("pagenum") != null) ? getIntParam("pagenum") : 1;
			Integer pageSize = (getParam("pagesize") != null) ? Math.min(getIntParam("pagesize"), 500) : 10;
			BatchUpoadBrowserFilter filter = new BatchUpoadBrowserFilter(getParam("type"), getParam("rownum"));
			List<RecordPreview> preview = processor.getRecordPreviews((pageNum - 1)*pageSize, pageSize, filter);
			WebList list = new WebList();
			list.list = (List)preview;
			list.pageNum = pageNum;
			list.pageSize = pageSize;
			list.size = processor.getUploadPreview().getRecords(filter).size();
			return new WebModel(list).getModelAndView();
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}

	public ModelAndView status(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		MultithreadBatchUploadProcessor processor = getProcessor();
		try
		{
			if (assertParam("interrupt"))
				processor.context.interruptRequested = true;

			String	status = processor.eh.message + (processor.context.interruptRequested ? " (interrupting)":"");

			return new WebModel(Alert.custom("status", processor.isThreadRunning() ? "" : getCurrentPage(), status)).setTemplate("batch30/status").getModelAndView();
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}

	public ModelAndView browser_action(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		if (!assertParam("action")) 
			return null;

		BatchUploadProcessor processor = getProcessor();
		try 
		{
			UploadPreview up = processor.getUploadPreview();

			if (up == null)
				return new WebModel().getModelAndView();

			if (req.getParameter("action").equals("state"))
			{
				for (Object param : req.getParameterMap().keySet())
				{
					String ps = (String)param;
					if (!ps.startsWith("upload_action"))
						continue;
					String value = getParam(ps);
					Long index = Long.valueOf(ps.replaceAll("upload_action", ""));
					for (RecordPreview rp : up)
						if (rp.id == index)
							rp.action = PreviewUploadAction.valueOf(value);
				}
			} else
				if (req.getParameter("action").equals("batchstate"))
				{
					String operation = getParam("boperation");
					String type = getParam("btype");
					for (RecordPreview rp : up)
						if (rp.status == PreviewRecordStatus.valueOf(type))
							rp.action = PreviewUploadAction.valueOf(operation);
				}

			return new WebModel().getModelAndView();
		} catch (Exception e)
		{
			return new WebModel(Alert.Exception(e)).setTemplate("batch30/exception").setRenderMode("popup").getModelAndView();
		}
	}
}