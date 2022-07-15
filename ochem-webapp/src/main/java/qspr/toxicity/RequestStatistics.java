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

package qspr.toxicity;

import java.util.Calendar;

import org.hibernate.stat.Statistics;

import qspr.Globals;

public class RequestStatistics
{
	Statistics hibernateStatistics;
	public long time;
	public long queryCount;
	public long insertCount;
	public long updateCount;
	public long deleteCount;
	
	public RequestStatistics()
	{
		hibernateStatistics = Globals.sessionFactory.getStatistics();
	}
	
	public void start()
	{
		time = Calendar.getInstance().getTimeInMillis();
		queryCount = hibernateStatistics.getQueryExecutionCount();
		insertCount = hibernateStatistics.getEntityInsertCount();
		updateCount = hibernateStatistics.getEntityUpdateCount();
		deleteCount = hibernateStatistics.getEntityDeleteCount();
	}
	
	public void stop()
	{
		time = Calendar.getInstance().getTimeInMillis() - time;
		queryCount = hibernateStatistics.getQueryExecutionCount() - queryCount;
		insertCount = hibernateStatistics.getEntityInsertCount() - insertCount;
		updateCount = hibernateStatistics.getEntityUpdateCount() - updateCount;
		deleteCount = hibernateStatistics.getEntityDeleteCount() - deleteCount;
	}
	
	public String current()
	{
		stop();
		String st = this.toString();
		start();
		
		return st;
	}
	
	public String toString()
	{
		return time+" ms., CRUD="+insertCount+"+"+queryCount+"+"+updateCount+"+"+deleteCount;
	}
}
