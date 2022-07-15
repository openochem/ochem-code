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

package qspr.metaserver;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import qspr.metaserver.protocol.ServerInfo;
import qspr.metaserver.protocol.Task;

/**
 * Class to cache info about recently active servers/client
 * 
 * @author midnighter
 */
public class OnlinePeer
{
	public long lastActive;
	public boolean isClient = false;
	public Task currentTask;
	public String ipAddress;
	public String status;
	
	/**
	 * A command to send to the peer next time it connects
	 */
	public String commandToSend;
	
	public ServerInfo serverInfo;
	public long conflictTime = 0;
	
	// This server was manually disabled
	public boolean disabled = false;
	
	// Some tasks from this server have been manually disabled
	public Set<String> disabledTasks = new HashSet<String>();
	
	public boolean restartRequested = false;
	public boolean logsRequested = false;
	
	public void disableTask(String taskType)
	{
		if (serverInfo.supportedTaskTypes.contains(taskType))
			disabledTasks.add(taskType);
	}
	
	public void enableTask(String taskType)
	{
		disabledTasks.remove(taskType);
	}
	
	public long getPing() // in milliseconds
	{
		long curTime = Calendar.getInstance().getTimeInMillis();
		return curTime - lastActive;
	}
	
	public boolean notResponding()
	{
		return getPing() > MetaServer.METASERVER_TASK_NOT_RESPONING;
	}
	
	public List<String> getAvailableTaskTypes()
	{
		List<String> availableTaskTypes = new ArrayList<String>();
		for (String availableTask : serverInfo.supportedTaskTypes)
			if (!disabledTasks.contains(availableTask))
				availableTaskTypes.add(availableTask);
		return availableTaskTypes;
	}
	
	
}
