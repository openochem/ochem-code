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

package qspr.workflow.utils;

import java.io.FileOutputStream;
import java.io.IOException;

import qspr.workflow.datatypes.DataTable;

public class RUtils 
{
	public static void datatableToRScript(DataTable dt, String filename) throws IOException
	{
		String st = RUtils.dataTableToMatrix(dt);
		FileOutputStream fos = new FileOutputStream(filename);
		fos.write(st.getBytes());
		fos.flush();
		fos.close();
	}
	
	public static String dataTableToMatrix(DataTable dt)
	{
		// copied from RServer. Maybe keep this method in one place
		StringBuffer buffer = new StringBuffer();
		buffer.append("datatable = matrix(c(");
		dt.reset();
		while (dt.nextRow())
		{
			for (int col = 0; col < dt.getColumnsSize(); col++)
				buffer.append(dt.getValue(col) + ", ");
		}
		buffer.delete(buffer.length()-2, buffer.length()-1);
		buffer.append("), "+dt.getRowsSize()+", "+dt.getColumnsSize()+", TRUE)");
		return buffer.toString();
	}
}
