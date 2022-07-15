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

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.frontend.WebModel;
import qspr.modelling.configurators.BasicModelConfigurator;
import qspr.util.Operation;
import qspr.util.StatusTracker;
import qspr.util.WrapperThread;

import com.eadmet.business.ModelUploadQuery;
import com.eadmet.business.ModelUploadResult;
import com.eadmet.business.ModelUploadService;


@Controller
public class ModelUploadController extends ControllerWrapper 
{
	//	private static transient final Logger logger = Logger.getLogger(ModelUploadController.class);

	public ModelUploadController()
	{
		sessionRequired = true;
	}

	public ModelAndView show(HttpServletRequest req, HttpServletResponse res) 
	{
		return new WebModel().setTemplate("modeller/model-upload").getModelAndView();
	}

	@SuppressWarnings("unchecked")
	public ModelAndView submit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		final ModelUploadQuery q = new ModelUploadQuery();
		q.trainingSetId = getLongParam("trainingsetid");
		q.validationSetIds = ModelController.getValidationSetIDs(request);
		q.unitMap = new HashMap<String, Long>();
		for (Object e : request.getParameterMap().entrySet()) 
		{
			Map.Entry<String, String[]> me = (Map.Entry<String, String[]>)e;
			if (me.getKey().matches("unit\\d+"))
			{
				System.out.println(me.getKey());
				System.out.println(Arrays.asList(me.getValue()));
				if (me.getValue().length > 0)
					if (me.getValue()[0].length() > 0)
						q.unitMap.put(me.getKey(), Long.valueOf(me.getValue()[0]));
			}
		}
		q.file = Globals.getUploadedFile("file");
		q.description = getParam("description");


		ModelUploadWrapperThread t = new ModelUploadWrapperThread()
		{
			@Override
			public void wrapped() throws Exception 
			{
				setStatus("Uploading model");
				ModelUploadService s = new ModelUploadService();
				OperationStatusChanger sc = new OperationStatusChanger();
				sc.o = operation;
				s.status.addListener(sc);
				result = s.uploadModel(q);
				setStatus("Finished");
			}
		};
		t.newOperation().start();
		Globals.setSessionAttribute(SessionVariable.MODEL_UPLOAD_RESULT, t);
		t.operation.successURL = "modelupload/results.do?operation-id=" + t.operation.operationId;
		return redirect("longoperations/operationWaitingScreen.do?operation-id=" + t.operation.operationId);
	}

	public ModelAndView results(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ModelUploadWrapperThread t = (ModelUploadWrapperThread)Globals.getSessionAttribute(SessionVariable.MODEL_UPLOAD_RESULT);
		for(int i=0; i<1000 && t.isAlive(); i++)
			Thread.sleep(100); // we need to wait until model is ready
		ModelUploadResult result = t.result;
		Globals.setSessionAttribute(SessionVariable.MODEL_UPLOAD_RESULT, null);
		if (result == null)
			throw t.exception;
		else
			if (result.errors.size() == 0)
			{
				BasicModelConfigurator configurator = new BasicModelConfigurator();
				configurator.setModel(result.model);

				Globals.setSessionAttribute(SessionVariable.MODEL_CONFIGURATOR, configurator);
				Globals.setSessionAttribute(SessionVariable.MODEL, result.model);

				return new WebModel(result.model).setTemplate("modeller/savemodel").addParam("page", "save").getModelAndView();
			}
			else
				return new WebModel(result).setTemplate("modeller/model-upload-error").getModelAndView();
	}
}

abstract class ModelUploadWrapperThread extends WrapperThread
{
	ModelUploadResult result;
}

class OperationStatusChanger extends StatusTracker
{
	Operation o;

	@Override
	public void set(String status)
	{
		super.set(status);
		if (!"Finished".equals(status))
			o.setStatus(status);
	}
}




