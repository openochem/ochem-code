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
import java.io.PrintWriter;
import java.util.List;

/**
 * A writer to export OCHEM data into the comma-separated format (CSV)
 * @author midnighter
 *
 */
public class CSVExportWriter extends ExportWriter
{
	public PrintWriter pw;

	@SuppressWarnings("rawtypes")
	@Override
	public void writeRow(List row)
	{
		StringBuffer line = new StringBuffer();
		for (int i = 0; i < row.size(); i++)
		{
			String val = row.get(i) == null ? "-" : "" + row.get(i);
			line.append("\"" + val + "\"");
			if (i < row.size() - 1)
				line.append(",");
		}

		pw.println(line);
	}

	@Override
	public void initialize()
	{
		pw = new PrintWriter(os);
		writeRow(columns);
	}


	@Override
	public void flush() throws IOException
	{
		pw.flush();
		super.flush();
	}

	@Override
	public void writeSupplementaryData(String key, Object value) 
	{
		// Unsupported for CSV yet
	}

	@Override
	public String getFileExtension() {
		return "csv";
	}

	@Override
	public String getHttpContentType() {
		return "application/csv";
	}

}
