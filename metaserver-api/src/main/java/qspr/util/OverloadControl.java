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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * A class that encapsulates the "waiting" logic to prevent overloading of misc. resources
 * Later we can add "incremental waiting" logic here.
 * 
 * Inspired by "deep sleep" concept.
 * 
 * @author midnighter&novserj
 *
 */
public class OverloadControl
{
	private static final Logger logger = LogManager.getLogger(OverloadControl.class);
	private long currentTimeout = 1000;
	private long timeout = 1000;
	private long maxTimeout = 10000;
	private long numMaxTimeout = 0;
	private long numTimeout = 0;
	
	public String resource;
	
	public void relax(Throwable e) throws InterruptedException
	{
		long waitTime = Math.round(currentTimeout * 2 * Math.random());
		logger.warn("Waiting " + waitTime + "ms. not to overload " + resource, e);
		Thread.sleep(waitTime);
		
		currentTimeout *= 1.5;
		if (currentTimeout > maxTimeout)
		{
			currentTimeout = maxTimeout;
			numMaxTimeout++;
		}
		numTimeout++;
	}
	
	public void reset()
	{
		currentTimeout = timeout;
	}
	
	public long getNumTimeout()
	{
		return numTimeout;
	}
	
	public long getNumMaxTimeout()
	{
		return numMaxTimeout;
	}
	
	public long getCurrentTimeout()
	{
		return currentTimeout;
	}
	
	public OverloadControl(String resource, long timeout, long maxTimeout)
	{
		this.resource = resource;
		this.timeout = timeout;
		this.currentTimeout = this.timeout;
		this.maxTimeout = maxTimeout;
	}
}
