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

import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.Logger;


/**
 * A simple incapsulation of the operation status.
 * Can simultaneously update status of a linked listener and print it to the logger.
 * 
 * @author midnighter
 *
 */
public class StatusTracker
{
	private String status;
	
	/**
	 * A linked status object updated simultaneously
	 */
	private List<StatusTracker> listeners;
	
	private Logger logger;
	private String prefix;
	
	public void set(String newStatus)
	{
		if (prefix != null)
			newStatus = prefix + newStatus;
		if (listeners != null)
			for (StatusTracker listener : listeners)
				listener.set(newStatus);
		if (logger != null && newStatus != null && !newStatus.equals(status))
			logger.info(newStatus);
		
		this.status = newStatus;
	}
	
	public String get()
	{
		return status;
	}
	
	public StatusTracker()
	{
		
	}
	
	public StatusTracker(Logger logger)
	{
		this.logger = logger;
	}
	
	public StatusTracker addListener(StatusTracker listener)
	{
		if (listeners == null)
			listeners = new ArrayList<StatusTracker>();
		listeners.add(listener);
		
		return this;
	}
	
	public StatusTracker setPrefix(String prefix)
	{
		this.prefix = prefix;
		return this;
	}
}
