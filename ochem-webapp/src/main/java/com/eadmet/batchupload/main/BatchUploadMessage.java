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

package com.eadmet.batchupload.main;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import qspr.Globals;

import com.eadmet.utils.MAILERConstants;
import com.eadmet.utils.OCHEMUtils;

@XmlRootElement
public class BatchUploadMessage
{
	@XmlElement
	public BatchUploadMessageType type;

	@XmlElement
	public String message;

	@XmlElement
	public String trace;

	public static BatchUploadMessage newError(String message)
	{
		BatchUploadMessage m = new BatchUploadMessage();
		m.type = BatchUploadMessageType.error;
		m.message = message;
		return m;
	}

	public static BatchUploadMessage newWarning(String message)
	{
		BatchUploadMessage m = new BatchUploadMessage();
		m.type = BatchUploadMessageType.warning;
		m.message = message;
		return m;
	}
	public static BatchUploadMessage newNotice(String message)
	{
		BatchUploadMessage m = new BatchUploadMessage();
		m.type = BatchUploadMessageType.notice;
		m.message = message;
		return m;
	}

	public static BatchUploadMessage newError(Throwable e)
	{
		BatchUploadMessage m = new BatchUploadMessage();
		m.type = BatchUploadMessageType.error;
		if (e instanceof NullPointerException)
			m.message = "Null pointer exception";
		else
			m.message = e.getMessage();
		if(Globals.userSession().user.login.equals(MAILERConstants.ADMIN))
			m.trace = OCHEMUtils.exceptionToString(e);
		return m;
	}

	public enum BatchUploadMessageType
	{
		error, warning, notice
	}

	public String toString()
	{
		return type+":"+message+"\n"+trace;
	}

}
