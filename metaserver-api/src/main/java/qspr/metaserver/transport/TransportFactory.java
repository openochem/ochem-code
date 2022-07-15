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

package qspr.metaserver.transport;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * Allows to configure the transport class used by default and to transparently inject non-standard implementations.
 * 
 * @author midnighter
 *
 */
public class TransportFactory 
{
	/**
	 * The default transport class for the whole JVM
	 */
	@SuppressWarnings("rawtypes")
	public static Class defaultTransportClass = CSTransport.class;

	/**
	 * The default transport for the current thread
	 */
	private static ThreadLocal<Transport> threadTransport = new ThreadLocal<Transport>();

	/**
	 * An interval (in seconds) which a calculation server should wait between task status requests.
	 * Not yet all calculation servers are using this value at the moment. 
	 */
	public static int waitingTimeBetweenRequests = 10 * 1000;

	public static Transport create()
	{
		try
		{
			if (threadTransport.get() != null)
			{
				return (Transport) threadTransport.get();
			}
			else
			{
				return (Transport) defaultTransportClass.newInstance();
			}
		}
		catch (Exception e)
		{
			throw new RuntimeException("Cannot create calculation server transport");
		}
	}

	public static void setThreadTransport(Transport t)
	{
		logger.debug("Setting thread-bounded Transport to " + t);
		threadTransport.set(t);
	}

	public static void clearThreadTransport()
	{
		logger.debug("Clearing thread-bounded transport");
		threadTransport.remove();
	}

	public static Transport getThreadTransport()
	{
		return threadTransport.get();
	}


	private static final Logger logger = LogManager.getLogger(TransportFactory.class);
}
