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
import java.util.zip.GZIPOutputStream;

/**
 * Allows to export OCHEM data into SDF format. 
 * The columns are exported as supplementary data.
 * @author midnighter
 *
 */
public class SDFExportWriter extends ExportWriter
{
	private GZIPOutputStream zipStream;
	private PrintWriter pw;

	@Override
	public boolean skipDescriptor(Object desc){
		if( desc == null)return true;
		if( desc instanceof Double && ((Double) desc == 0D))return true;
		return false;
	}


	@Override
	public void initialize() throws IOException
	{
		zipStream = new GZIPOutputStream(os);
		pw = new PrintWriter(zipStream);
	}

	@SuppressWarnings("rawtypes")
	@Override
	public void writeRow(List row)
	{
		pw.println(cleanupSDF(currentMolecule.sdf));
		for (int i = 0; i < row.size(); i++)
		{
			Object val = row.get(i);
			if (val != null && !("" + val).isEmpty())
			{
				pw.println("> <"+columns.get(i)+">");
				pw.println("" + val);
				pw.println();
			}
		}
		pw.println("$$$$");
	}

	private String cleanupSDF(String sdf)
	{
		return sdf.replaceAll("M  END[\\s\\S]*", "M  END");
	}

	@Override
	public void flush() throws IOException
	{
		pw.flush();
		zipStream.flush();
		pw.close();
		zipStream.close();

		super.flush();
	}

	@Override
	public void writeSupplementaryData(String key, Object value) 
	{
		// Unsupported for SDF yet
	}

	@Override
	public String getFileExtension() {
		return "sdf.gz";
	}

	@Override
	public String getHttpContentType() {
		return "application/x-gzip";
	}

}
