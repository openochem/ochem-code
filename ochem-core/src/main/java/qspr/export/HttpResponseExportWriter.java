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

package qspr.export;

import java.io.IOException;
import java.util.List;

import javax.servlet.http.HttpServletResponse;

/**
 * A writer to export OCHEM data directly into HTTP stream using a downstream format-specific writer
 * @author novserj
 *
 */


public class HttpResponseExportWriter extends ExportWriter
{
	protected HttpServletResponse response;
	protected ExportWriter underlyingWriter;

	public void writeSupplementaryData(String key, Object value)
	{
		underlyingWriter.writeSupplementaryData(key, value);
	}

	@SuppressWarnings("rawtypes")
	public void writeRow(List row)
	{
		underlyingWriter.writeRow(row);
	}

	public void initialize() throws IOException
	{
		underlyingWriter.os = response.getOutputStream();
		response.setContentType(underlyingWriter.getHttpContentType());
		response.setHeader("Content-Disposition", "attachment; filename="+underlyingWriter.fileName+"."+underlyingWriter.getFileExtension());
	}

	@Override
	public void setFileName(String fileName)
	{
		underlyingWriter.setFileName(fileName);
	}

	@Override
	public void write() throws Exception
	{
		initialize();
		underlyingWriter.write();
	}

	public HttpResponseExportWriter(HttpServletResponse response, ExportWriter underlyingWriter)
	{
		this.underlyingWriter = underlyingWriter;
		this.response = response;
	}

	@Override
	public String getFileExtension() {
		return null;
	}

	@Override
	public String getHttpContentType() {
		return null;
	}

}
