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
import java.lang.reflect.Method;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.springframework.web.servlet.ModelAndView;

import qspr.entities.Alert;
import qspr.frontend.WebModel;
import qspr.tests.PeriodicTestRunner;
import qspr.util.Operation;

public class StaticController extends ControllerWrapper {
	
	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) 
	throws Exception 
	{
		String[] pieces = request.getRequestURI().split("/");
		pieces = pieces[pieces.length-1].split("\\.");
		
		// Is there a method for this request here?
		try
		{
			Method method = this.getClass().getMethod(pieces[0], HttpServletRequest.class, HttpServletResponse.class);
			return (ModelAndView) method.invoke(this, request, response);
		}
		catch (NoSuchMethodException e)
		{
			WebModel wm = new WebModel();
			wm.templateName = pieces[0];
			return wm.getModelAndView();
		}
	}
	
	public ModelAndView fileDownloadScreen(HttpServletRequest request, HttpServletResponse response) 
	throws Exception 
	{
		String name = getParam("filename");
		return new WebModel(new Alert(name)).setTemplate("file-download").getModelAndView();
	}
	
	public ModelAndView shutdown(HttpServletRequest request, HttpServletResponse response) 	throws Exception 
	{
		System.exit(0);
		return null;
	}
	
	/**
	 * Force running server tests for debugging/coverage purposes
	 */
	public ModelAndView runServerTests(HttpServletRequest request, HttpServletResponse response) 
			throws Exception 
	{
		PeriodicTestRunner runner = new PeriodicTestRunner();
		PeriodicTestRunner.loader = this.getClass().getClassLoader();
		runner.filter.disableScheduler = true;
		runner.listener.saveToDB = false;
		runner.run();
		return null;
	}
	
	// Universal method to retrieve the status of a long running operation
	public ModelAndView operationStatus(HttpServletRequest request, HttpServletResponse response) 
	throws Exception 
	{
		Long operationId = getLongParam("operation-id");
		Operation operation = Operation.getOperation(operationId);
		String status;
		if (operation != null)
			status = operation.getStatus();
		else
			status = "Error: Unknown operation. Probably, the server has been restarted during the operation.\nPlease, retry - it should not happen again.";
		
		return new WebModel(new Alert(status)).getModelAndView();
	}
}
