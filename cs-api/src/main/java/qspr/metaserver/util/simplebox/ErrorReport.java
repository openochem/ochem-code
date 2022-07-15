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

package qspr.metaserver.util.simplebox;

import java.io.IOException;
import java.io.Writer;
import java.util.HashMap;
import java.util.LinkedHashMap;

import org.json.simple.JSONValue;

public class ErrorReport extends AbstractComponent {

	/** The HTTP status that accompanied the Error Report */
	private int httpStatus;
	/** The peer that threw the exception or reported an exceptional event */
	private String actor;
	/** Brief explanatory message */
	private String message;
	/** Technical Details */
	private String details;
	/** Error Cause Identification Code */
	private String errorCode;
	/** Trace... */
	private ErrorReport errorCause;

	public ErrorReport() {
	}

	public ErrorReport(int httpStatus, String actor, String message,
			String details, String errorCode) {
		this.httpStatus = httpStatus;
		this.actor = actor;
		this.message = message;
		this.details = details;
		this.errorCode = errorCode;
	}

	public String getActor() {
		return actor;
	}

	public void setActor(String actor) {
		this.actor = actor;
	}

	public String getDetails() {
		return details;
	}

	public void setDetails(String details) {
		this.details = details;
	}

	public ErrorReport getErrorCause() {
		return errorCause;
	}

	public void setErrorCause(ErrorReport errorCause) {
		this.errorCause = errorCause;
	}

	public int getHttpStatus() {
		return httpStatus;
	}

	public void setHttpStatus(int errorCode) {
		this.httpStatus = errorCode;
	}

	public String getMessage() {
		return message;
	}

	public void setMessage(String message) {
		this.message = message;
	}

	public String getErrorCode() {
		return errorCode;
	}

	public void setErrorCode(String errorCode) {
		this.errorCode = errorCode;
	}

	@Override
	public void writeJSONString(Writer writer) throws IOException {

		HashMap<String, Object> obj = new LinkedHashMap<String, Object>();
		obj.put("object_type", "ErrorReport");
		if (getErrorCode() != null) {
			obj.put("code", getErrorCode());
		}
		if (message != null) {
			obj.put("message", message);
		}

		if (getDetails() != null) {
			obj.put("debug", getDetails());
		}

		if (getHttpStatus() > 0) {
			obj.put("httpcode", Integer.valueOf(getHttpStatus()));
		}

		if (getActor() != null) {
			obj.put("toBlame", getActor());
		}

		if (getErrorCause() != null) {
			obj.put("trace", getErrorCause());
		}

		if (getMeta() != null) {
			obj.put("metadata", getMeta());
		}
		JSONValue.writeJSONString(obj, writer);
	}

}
