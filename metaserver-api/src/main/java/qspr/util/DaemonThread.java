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

import org.apache.logging.log4j.Logger;

public abstract class DaemonThread extends Thread
{
	abstract public void wrapped() throws Exception;

	protected Logger log;
	
	public DaemonThread(Logger logger)
	{
		super();
		log = logger;
		setName(this.getClass().getSimpleName());
		setDaemon(true);
	}

	@Override
	public void run()
	{
		try
		{
			while (true)
			{
				wrapped();
			}
		} catch (InterruptedException e)
		{
			log.info("The daemon thread (" + this.getClass().getSimpleName() + ") has been interrupted");
		} catch (Exception e)
		{
			log.info("SEVERE! The daemon thread (" + this.getClass().getSimpleName() + ") has exited with an exception: " + e.getMessage());
		}
	}
}
