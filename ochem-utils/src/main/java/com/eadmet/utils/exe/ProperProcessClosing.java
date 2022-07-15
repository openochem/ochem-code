/**
 * 
 */
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

package com.eadmet.utils.exe;

import java.io.Closeable;
import java.io.IOException;

/**
 * @author robert
 * 
 * according to java bugs
 * http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=4523660
 * http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=4801027
 * http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=6462165
 * 
 * I put this stuff everywhere a process is called, that the streams will be closed 
 * properly and "Cannot run program "[some*bash*program]": java.io.IOException: error=24, Too many open files"
 * will occur anymore. Suggested by
 * http://stuffthathappens.com/blog/2007/11/28/crash-boom-too-many-open-files/
 *  
 */
public class ProperProcessClosing {

	public static void closeProcess(Process proc) {
		if (proc != null) {
	        close(proc.getOutputStream());
	        close(proc.getInputStream());
	        close(proc.getErrorStream());
	        proc.destroy();
	      }
	}
	
	private static void close(Closeable c) {
		if (c != null) {
	      try {
	        c.close();
	      } catch (IOException e) {
	        // ignored
	      }
	    }
	}
}