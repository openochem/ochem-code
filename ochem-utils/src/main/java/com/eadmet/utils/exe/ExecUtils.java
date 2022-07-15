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

import java.io.IOException;
import java.io.StringWriter;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Miscellaneous utilities related to shell commands
 * @author midnighter
 *
 */
public class ExecUtils
{
	/**
	 * Gets the output of a particular process as a string
	 */
	public static String getProcessOutput(Process p) throws InterruptedException
	{
		final StringWriter sw = new StringWriter();
		OutReader or = new OutReader(p.getInputStream(), Thread.currentThread()){
			public void onLineRead(String line) throws Exception {
				sw.append(line + "\n");
			};
		};
		
		OutReader er = new OutReader(p.getErrorStream(), Thread.currentThread()){
			public void onLineRead(String line) throws Exception {
				sw.append(line + "\n");
			};
		};
		or.start();
		er.start();
		p.waitFor();
		or.join();
		er.join();
		
		return sw.toString();
	}
	
	/**
	 * Gets the output of a particular shell command as a string
	 * @throws IOException 
	 */
	public static String getShellOutput(String shellCommand) throws InterruptedException, IOException
	{
		logger.info("Running command " + shellCommand);
		String[] commands = { "/bin/bash", "-c", shellCommand };
		Process p = Runtime.getRuntime().exec(commands, null);
		return getProcessOutput(p);
	}
	
	private static final Logger logger = LoggerFactory.getLogger(ExecUtils.class);
	
	public static void main(String[] args) throws IOException, InterruptedException
	{
		String[] commands = { "/bin/bash", "-c", "/Applications/Calculator.app/Contents/MacOS/Calculator 2>&1 > /dev/null &" };
		Process p = Runtime.getRuntime().exec(commands, null);
		Thread.sleep(10000);
		p.destroy();
	}
}
