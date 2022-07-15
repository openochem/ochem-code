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
import java.util.ArrayList;
import java.util.List;

import qspr.util.RWriter;

/**
 * Allows to export OCHEM data into R format
 * @author midnighter
 *
 */

public class RExportWriter extends ExportWriter
{
	RWriter writer;
	List<String> validColumnIds = new ArrayList<String>();

	@SuppressWarnings("rawtypes")
	@Override
	public void writeRow(List row)
	{
		for (int i = 0; i < row.size(); i++)
			writer.addValue(validColumnIds.get(i), row.get(i));
	}

	@Override
	public void initialize() throws IOException
	{
		writer = new RWriter(os);
		writer.variablePrefix = "model$";

		for (String column : columns)
			validColumnIds.add(column.replaceAll("[^a-zA-Z0-9]+", ".").replaceAll("\\.+$", ""));
	}

	@Override
	public void flush() throws IOException
	{
		writer.write();
		super.flush();
	}

	@Override
	public void writeSupplementaryData(String key, Object value) 
	{
		// Unsupported for R yet
	}

	@Override
	public String getFileExtension() {
		return "R";
	}

	@Override
	public String getHttpContentType() {
		return "application/R";
	}
}
