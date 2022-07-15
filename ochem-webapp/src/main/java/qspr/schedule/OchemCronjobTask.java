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

package qspr.schedule;

import java.util.Calendar;

import org.apache.logging.log4j.LogManager;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.util.Logger;

import com.eadmet.utils.mailer.Mailer;

/**
 * 
 * Abstract superclass for all Quartz-based OCHEM cron jobs
 * @author midnighter
 *
 */

public abstract class OchemCronjobTask implements Logger
{
	protected static transient final org.apache.logging.log4j.Logger logger = LogManager.getLogger(OchemCronjobTask.class);

	public OchemCronjobTask()
	{

	}

	/**
	 * Defines whether this cron job should run at this particular installation.
	 * Can be overriden or complemented by subclasses
	 */
	protected boolean shouldRun()
	{
		if(OCHEMConfiguration.mirror) return false; 

		if (this.getClass().getAnnotation(DatabaseMaintenanceJob.class) != null)
			return !OCHEMConfiguration.testing;

		if (this.getClass().getAnnotation(GlobalMaintenanceJob.class) != null)
			return !OCHEMConfiguration.testing;

		return true;
	}

	public void executeInternal()
	{
		if (shouldRun())
		{
			try
			{
				log("Task started");
				long time = Calendar.getInstance().getTimeInMillis();
				executeTask();
				log("Task finished in "+(Calendar.getInstance().getTimeInMillis() - time)+"ms");
			} catch (Exception e)
			{
				log("Task failed");
				e.printStackTrace();
				Globals.rollbackAllTransactions();
				Mailer.notifyDevelopers(e, this.getClass().getName());
			}
		} 
		//else
		//	log("Skipping");
	}

	public void log(String st)
	{
		String className = this.getClass().getName();
		className = className.substring(className.lastIndexOf(".")+1);
		logger.info("["+Globals.now()+"]["+className+"] " + st);
	}

	public abstract void executeTask() throws Exception;
}
