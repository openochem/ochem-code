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

import java.util.Calendar;

import org.apache.logging.log4j.Logger;

public class SmartLogger
{
	private Logger logger;
	private long time;
	private long delta = 5000;
	private StatusTracker statusTracker;
	
	public SmartLogger(Logger logger) {
		this.logger = logger;
		time = Calendar.getInstance().getTimeInMillis();
	}
	
	public SmartLogger(Logger logger, long intervalInMilliseconds) {
		this.logger = logger;
		this.delta = intervalInMilliseconds;
		time = Calendar.getInstance().getTimeInMillis();
	}
	
	public SmartLogger setStatusTracker(StatusTracker tracker) {
		this.statusTracker = tracker;
		return this;
	}
	
	public void log(String msg) {
		if (itsTime())
		{
			logger.info(msg);
			if (statusTracker != null)
				statusTracker.set(msg);
		}
	}
	
	/**
	 * Check if the specified interval passed since the last call
	 * Reset the timer if the interval has passed
	 */
	public boolean itsTime() {
		long newTime = Calendar.getInstance().getTimeInMillis();
		if (newTime - time > delta) {
			time = newTime;
			return true;
		}
		
		return false;
	}
}
