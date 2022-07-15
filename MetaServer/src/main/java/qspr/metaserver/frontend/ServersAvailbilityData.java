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

package qspr.metaserver.frontend;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import qspr.metaserver.MetaServer;
import qspr.metaserver.OnlinePeer;

/**
 * A summary of servers availability and task queue by task types and priorities
 * @author midnighter
 *
 */
public class ServersAvailbilityData
{
	public Map<String, TaskTypeAvailability> availabilityByTaskType = new HashMap<String, ServersAvailbilityData.TaskTypeAvailability>();
	
	/**
	 * All available task priorities
	 */
	public Set<Integer> priorities = new HashSet<Integer>();
	
	/**
	 * All available task types
	 */
	public Set<String> taskTypes = new HashSet<String>();
	
	/**
	 * Get count of (free) servers for a task type and priority
	 */
	public int getAvailability(String taskType, Integer priority, boolean free)
	{
		TaskTypeAvailability tta = availabilityByTaskType.get(taskType);
		if (tta == null)
			return 0;
		return free ? tta.getAvailability(priority).freeServers : tta.getAvailability(priority).totalServers;
	}
	
	/**
	 * Gather the servers availability data based on the list of online peers
	 * @param peers
	 */
	public void addPeers(Collection<OnlinePeer> peers)
	{
		for (OnlinePeer peer : peers)
		{
			if (peer.isClient)
				continue;
			if (peer.serverInfo == null)
				continue;
			
			if (peer.getPing() > MetaServer.METASERVER_TASK_NOT_RESPONING)
				continue;
			
			for (String taskType : peer.serverInfo.supportedTaskTypes)
				addTaskType(taskType, peer.serverInfo.minimumPriority, peer.currentTask == null);
		}
	}
	
	public void addTask(String taskType, Integer priority, boolean assigned)
	{
		if (priority == null)
			priority = -1000000;
		priorities.add(priority);
		if (assigned)
			getAvailability(taskType, priority).inCalculation++;
		else
			getAvailability(taskType, priority).queue++;
	}
	
	public TaskTypePriorityInfo getAvailability(String taskType, Integer priority)
	{
		TaskTypeAvailability tta = getAvailability(taskType);
		return tta.getAvailability(priority);
	}
	
	private TaskTypeAvailability getAvailability(String taskType)
	{
		TaskTypeAvailability tta = availabilityByTaskType.get(taskType);
		if (tta == null)
		{
			taskTypes.add(taskType);
			availabilityByTaskType.put(taskType, tta = new TaskTypeAvailability());
		}
		return tta;
	}
	
	public void addTaskType(String taskType, Integer priority, boolean free)
	{
		if (priority == null)
			priority = -1000000;
		priorities.add(priority);
		TaskTypeAvailability tta = getAvailability(taskType);
		TaskTypePriorityInfo availability = tta.getAvailability(priority);
		if (free)
			availability.freeServers++;
		availability.totalServers++;
	}
	
	public static class TaskTypeAvailability
	{
		/**
		 * Maps the priority to the two {free, total} number of available servers
		 */
		public Map<Integer, TaskTypePriorityInfo> availabilityByPriority = new HashMap<Integer, TaskTypePriorityInfo>();
		
		public TaskTypePriorityInfo getAvailability(Integer priority)
		{
			TaskTypePriorityInfo res = availabilityByPriority.get(priority);
			if (res == null)
			{
				res = new TaskTypePriorityInfo();
				availabilityByPriority.put(priority, res);
			}
			
			return res;
		}
	}
	
	/**
	 * Information about free+total servers and queued tasks for a particular task type and priority
	 * @author midnighter
	 *
	 */
	public static class TaskTypePriorityInfo {
		public int freeServers;
		public int totalServers;
		public int queue;
		public int inCalculation;
	}
}
