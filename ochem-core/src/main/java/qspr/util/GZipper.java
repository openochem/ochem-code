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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class GZipper
{

	public static OutputStream fileToOutputStream(String fileName)
	{
		GZIPOutputStream gzipOS = null;

		if (!fileName.endsWith(".gz"))
		{
			fileName += ".gz";
		}

		try
		{
			FileOutputStream fOS = new FileOutputStream(fileName);
			gzipOS = new GZIPOutputStream(fOS);
		} 
		catch (IOException e)
		{
			e.printStackTrace();
		}

		return gzipOS;
	}
	
	public static InputStream fileToInputStream(String fileName)
	{
		FileInputStream fIS = null;
		GZIPInputStream gzipIS = null;
		
		try 
		{
			if ( ! fileName.endsWith(".gz"))
			{
				File f = new File(fileName);
				if (f.exists())
					return new FileInputStream(fileName);
				else
					fileName += ".gz";
			}

			fIS = new FileInputStream(fileName);
			gzipIS = new GZIPInputStream(fIS);
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} 

		return gzipIS;
	}
}
