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

package qspr.util;

import javax.servlet.http.HttpServletRequest;

import org.springframework.web.multipart.MultipartHttpServletRequest;

import qspr.ThreadScope;

/**
 * A convenience class to parse HTTP requests
 * @author midnighter
 *
 */
public class RequestParser
{
	HttpServletRequest request;

	public RequestParser(HttpServletRequest request)
	{
		this.request = request;
	}

	public RequestParser()
	{
		MultipartHttpServletRequest mpRequest = ThreadScope.get().localMpRequest;
		if (mpRequest == null)
			this.request = ThreadScope.get().localRequest;
		else
			this.request = mpRequest;
	}

	public boolean assertParam(String name)
	{
		return request.getParameter(name) != null && !request.getParameter(name).equals("") && !request.getParameter(name).equals("undefined");
	}

	public Long getLongParam(String name)
	{
		if (!assertParam(name))
			return null;
		String param = request.getParameter(name);
		if (param != null)
			return Long.valueOf(param);
		else
			return null;

	}

	public Integer getIntParam(String name)
	{
		if (!assertParam(name))
			return null;
		return Integer.valueOf(request.getParameter(name));
	}

	public String getParam(String name)
	{
		return request.getParameter(name);
	}
}
