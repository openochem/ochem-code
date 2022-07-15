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

package qspr.services;

import com.eadmet.utils.OCHEMUtils;

public class CreateEntityResponse
{
	private Long entityId;
	private String status;
	private String message;
	
	public Long getEntityId()
	{
		return entityId;
	}

	public void setEntityId(Long entityId)
	{
		this.entityId = entityId;
	}

	public String getStatus()
	{
		return status;
	}

	public void setStatus(String status)
	{
		this.status = status;
	}

	public String getMessage()
	{
		return message;
	}

	public void setMessage(String message)
	{
		this.message = message;
	}

	public CreateEntityResponse(long id)
	{
		this.entityId = id;
		this.status = "ok";
	}
	
	public static CreateEntityResponse ok(String message)
	{
		CreateEntityResponse response = new CreateEntityResponse(0);
		response.status = "ok";
		response.message = message;
		return response;
	}
	
	public static CreateEntityResponse status(String status, String message)
	{
		CreateEntityResponse response = new CreateEntityResponse(0);
		response.status = status;
		response.message = message;
		return response;
	}
	
	public CreateEntityResponse(Exception e)
	{
		this.status = "error";
		this.message = OCHEMUtils.exceptionToString(e);
	}
	
	public CreateEntityResponse(String errorMessage)
	{
		this.status = "error";
		this.message = errorMessage;
	}
}
