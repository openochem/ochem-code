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

import qspr.entities.Alert;
import qspr.frontend.WebModel;
import qspr.util.Operation;

@Controller
public class LongOperationsController extends ControllerWrapper
{
	public void show(HttpServletRequest request, HttpServletResponse response)
	{
		
	}
	
	public ModelAndView operationWaitingScreen(HttpServletRequest request, HttpServletResponse response) 
	throws Exception 
	{
		return new WebModel().setTemplate("operation-wait").getModelAndView();
	}
	
	public ModelAndView operationSuccess(HttpServletRequest request, HttpServletResponse response) 
	throws Exception 
	{
		Long operationId = getLongParam("operation-id");
		Operation operation = Operation.getOperation(operationId);
		return redirect(operation.successURL);
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
	
	public ModelAndView listOperations(HttpServletRequest request, HttpServletResponse response)
	{
		WebModel wm = new WebModel();
		for (Operation operation : Operation.operations.values())
			wm.addObject(operation);
				
		return wm.setTemplate("operations").getModelAndView();
	}
	
	public ModelAndView cancelOperation(HttpServletRequest request, HttpServletResponse response) 
	throws Exception 
	{
		Long operationId = getLongParam("operation-id");
		Operation operation = Operation.getOperation(operationId);
		if (operation != null)
			operation.cancel();
		
		return new WebModel().getModelAndView();
	}
}
