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

package qspr.metaserver.util.simplebox.parsers;

import org.json.simple.JSONObject;

import qspr.metaserver.util.simplebox.ErrorReport;

public class ErrorReportParser extends AbstractSBParser<ErrorReport> {

	public ErrorReportParser(JSONObject json) {
		super(json);
	}

	@Override
	public ErrorReport parse() throws JSONParsingException {
		assertObjectTypePresent("ErrorReport");
		ErrorReport er = new ErrorReport();

		Object code = json.get("code");
		er.setErrorCode(code != null ? code.toString() : null);

		Object httpStatus = json.get("httpcode");
		if (httpStatus != null) {
			try {
				int httpStatusInt = Integer.parseInt(httpStatus.toString());
				er.setHttpStatus(httpStatusInt);
			} catch (final NumberFormatException nfe) {
				// forgive this error - do nothing
				// TODO Log malformed JSON
			}
		}

		Object message = json.get("message");
		er.setMessage(message != null ? message.toString() : null);

		Object debug = json.get("debug");
		er.setDetails(debug != null ? debug.toString() : null);

		Object actor = json.get("toBlame");
		er.setActor(actor != null ? actor.toString() : null);

		Object trace = json.get("trace");
		if (trace != null) {
			JSONObject json_trace = (JSONObject) trace;
			ErrorReportParser traceParser = new ErrorReportParser(json_trace);
			ErrorReport mytrace = traceParser.parse();
			er.setErrorCause(mytrace);
		}

		attachMetadata(er);

		return er;
	}
}
