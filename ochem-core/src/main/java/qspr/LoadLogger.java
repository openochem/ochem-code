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

package qspr;

import java.io.File;

import com.eadmet.utils.exe.ProperProcessClosing;

public class LoadLogger 
{
	public static void log()
	{
		File f = new File(OCHEMConfiguration.customLogsPath+"/load.log");
		Process proc = null;
		try 
		{
			String top;
			
			if (System.getProperty("os.name").contains("Mac"))
				top = "top -l 1 | head -n 11";
			else
				top = "top -b -n 1 | head -n 6";
			
			proc = Runtime.getRuntime().exec(new String[]{"/bin/sh", "-c", top+"  >> "+f.getAbsolutePath()});
			proc.waitFor();
		} catch (Exception e) 
		{
			e.printStackTrace();
		} finally 
		{
			ProperProcessClosing.closeProcess(proc);
		}			
	}
}
