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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import qspr.metaserver.protocol.Task;

public class TaskQueue 
{
	public Map<String, List<Task>> tasksQueue = new HashMap<String, List<Task>>(); // "Workflow"->Task()

	public void clear() 
	{
		tasksQueue.clear();
	}

	public Task getBestTaskForPeer(OnlinePeer peer) 
	{
		Task taskToAssign = null;
		List<String> availableTasks = peer.getAvailableTaskTypes();
		for (String availableTask : availableTasks) 
		{
			List<Task> tasks = tasksQueue.get(availableTask);
			if (tasks == null)
				continue;
			for (Task task : tasks) // Possibly we can optimize to return the most demadning task, memory-wise
			{
				if ((peer.serverInfo.minimumPriority == null || peer.serverInfo.minimumPriority <= task.priority)
						&& (peer.serverInfo.availableMemory >= task.getMinRequiredMemory() * 0.9)) // Is the task important enough? Does the server have enough memory? Min required memory * 0.9 - to allow for incorrect calculation server configurations
				{
					if (taskToAssign == null || task.priority > taskToAssign.priority) // Do we have a more important task? (order of tasks in the list of supported tasks matters only if the priorities is equal) / Midnighter on Aug 8, 2011
						taskToAssign = task;
				}
			}
		}
		return taskToAssign;
	}

	private boolean tasksHaveSameMemoryRequirements(Task t1, Task t2) 
	{
		return  t1.getMinRequiredMemory() == t2.getMinRequiredMemory();
	}

	// We keep a list of most high priority tasks for each unique memory requirement
	// Therefore we put a task to a list if it has a unique memory requirement
	// or is higher priority than the task with the same memory requirement
	public boolean shouldPutTask(Task task) 
	{
		if (!Task.INIT.equals(task.status))
			return false;

		if (!tasksQueue.containsKey(task.taskType))
			return true;

		if (tasksQueue.get(task.taskType) == null)
			return true;

		if (tasksQueue.get(task.taskType).size() == 0)
			return true;

		List<Task> tasks = tasksQueue.get(task.taskType);
		for (int i = 0; i < tasks.size(); i++)
			if (tasksHaveSameMemoryRequirements(task, tasks.get(i)))
				return (task.hasHigherPriorityThan(tasks.get(i)));
		return true;
	}

	public void putTask(Task task) 
	{
		if (shouldPutTask(task)) 
		{
			List<Task> tasks = tasksQueue.get(task.taskType);

			if (tasks == null)
				tasksQueue.put(task.taskType, tasks = new ArrayList<Task>());

			for (int i = 0; i < tasks.size(); i++)
				if (tasksHaveSameMemoryRequirements(task, tasks.get(i))) 
				{
					if (task.hasHigherPriorityThan(tasks.get(i))) 
					{
						tasks.remove(i);
						tasks.add(task);
					}
					return;
				}
			tasks.add(task);
		}
	}

	public void removeTask(Task task) 
	{
		List<Task> tasks = tasksQueue.get(task.taskType);
		if (tasks != null)
			tasks.remove(task);
	}

	public Task getTaskById(int taskId) 
	{
		for (Entry<String, List<Task>> entry : tasksQueue.entrySet())
			if (entry.getValue() != null)
				for (Task task : entry.getValue())
					if (task.id.equals(taskId))
						return task;
		return null;
	}

	public boolean removeTaskById(int taskId) 
	{
		Task t = getTaskById(taskId);

		if (t != null)
			if (!t.hasMultipleReferences())
			{
				tasksQueue.get(t.taskType).remove(t);
				return true;
			}
		
		return false;
	}

}
