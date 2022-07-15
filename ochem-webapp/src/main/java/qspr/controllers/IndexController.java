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

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.VirtualHostConfiguration;

@Controller
public class IndexController extends ControllerWrapper 
{
	public IndexController()
	{
		useIframe = false;
		sessionRequired = false;
	}
	
	public ModelAndView index(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException
	{
		if (VirtualHostConfiguration.getHomePage() != null)
			return redirect(VirtualHostConfiguration.getHomePage() + "?render-mode=full");
		else
		{
			request.getRequestDispatcher("/index.html").forward(request, response);
			return null;
		}
	}
}
