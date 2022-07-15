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

import qspr.Globals;
import qspr.business.PendingTaskPeer;
import qspr.util.WrapperThread;

@GlobalMaintenanceJob
public class SynchronizePendingTasksTask extends OchemCronjobTask {

	@Override
	public void executeTask() throws Exception 
	{
		new WrapperThread() {
			@Override
			public void wrapped() throws Exception 
			{
				for (int i = 0; i < 100; i++)
				{
					boolean workDone = PendingTaskPeer.updateTaskStatuses(null, null, 100);
					Globals.restartAllTransactions(true);
					if (!workDone)
						break;
				}
			}
		}.run();
	}
	
	public static void main(String[] args) throws Exception
	{
		new SynchronizePendingTasksTask().executeTask();
	}
}
