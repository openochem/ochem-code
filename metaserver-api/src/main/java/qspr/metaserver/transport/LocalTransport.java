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

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.CalculationServer;
import qspr.metaserver.ServerPool;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.Task;

import com.eadmet.utils.MemoryUtils;

public class LocalTransport implements Transport
{
	private static final Logger logger = LogManager.getLogger(LocalTransport.class);
	public static Map<Integer, CalculationThread> runningTasks = new HashMap<Integer, CalculationThread>();
	public static Map<Integer, Task> queuedTasks = new HashMap<Integer, Task>();
	
	public boolean embeddedTransport = false;
	
	/**
	 * Allow queuing local tasks if the server is busy. Allows to calculate multiple tasks of the same type locally, e.g. for bagging
	 */
	public static boolean allowQueue = false;
	
	public synchronized Command executeCommand(Command command) throws UnsupportedOperationException 
	{
		Integer taskId;
		CalculationThread cThread;
		switch (command.id)
		{
			case Command.CL_SUBMIT_TASK:
				Task task = (Task) command.data;
				
				if (embeddedTransport && task.getMinRequiredMemory() > 0)
				{
					System.gc();
					long availableMemory = MemoryUtils.getCurrentMemoryFree();
					
					logger.info("Available memory in LocalTransport " + availableMemory + "MB required memory "+task.getMinRequiredMemory()+ "MB for "+task);
					
					availableMemory += (int) availableMemory * 0.2 <128?128:availableMemory * 0.2; // provides a tolerance, e.g. to run ASNNP inside of ASNN
					
					if (availableMemory < task.getMinRequiredMemory() && task.getMinRequiredMemory() > 1024) // for very small tasks we ignore estimation of memory; these tasks should run on one server!  
						throw new UnsupportedOperationException("Not enough free memory to run task" + task + " locally");
				}
				
				CalculationServer server = ServerPool.getInstance().getServer(task.taskType, true);
				
				// No free local servers. Try to get a busy one
				if (server == null && allowQueue)
					server = ServerPool.getInstance().getServer(task.taskType, false);
				
				if (server == null)
					throw new UnsupportedOperationException("Cannot run task type " + task.taskType + " locally");
				
				taskId = -Integer.valueOf((int)Math.round(Math.random()*1000000));
				
				synchronized (server) 
				{
					if (server.busy)
					{
						if (!allowQueue)
							throw new UnsupportedOperationException("Cannot run task type " + task.taskType + " locally");
						else
						{
							queuedTasks.put(taskId, task);
						}
					}
					server.busy = true;
				}
				
				task.id = taskId;
				task.calcServerId = ServerPool.getInstance().sid;
				task.status = Task.INIT;
				if (!queuedTasks.containsKey(taskId))
					runningTasks.put(taskId, new CalculationThread(server, task));
				return new Command(Command.MS_TASK_STATUS, task);
			case Command.CL_QUERY_TASK:
			case Command.CL_QUERY_TASK_STATUS:
				
				// Multiple statuses not yet supported
				if (command.data instanceof List)
					throw new UnsupportedOperationException();
				
				taskId = (Integer)command.data;
				
				synchronized (queuedTasks)
				{
					cThread = runningTasks.get(taskId);
					if (cThread == null)
						if (!queuedTasks.containsKey(taskId))
							throw new UnsupportedOperationException("Cannot request status of task " + taskId);
						else
							return new Command(Command.MS_TASK_STATUS, queuedTasks.get(taskId));
					else
					{
						if (cThread.task.isReady())
							runningTasks.remove(cThread);
						else
							cThread.task.setDetailedStatus(cThread.server.status);
						return new Command(Command.MS_TASK_STATUS, cThread.task);
					}
				}
			case Command.KILL_TASK:
				taskId = (Integer)command.data;
				cThread = runningTasks.get(taskId);
				if (cThread != null)
				{
					cThread.interrupt();
					runningTasks.remove(cThread);
				}
				else
					if (queuedTasks.containsKey(taskId))
						queuedTasks.remove(taskId);
					else
						throw new UnsupportedOperationException();
				
				return null;
			case Command.CL_GET_TASK_BY_MD5:
				if (embeddedTransport)
					throw new UnsupportedOperationException();
				return null;
		}
		
		throw new UnsupportedOperationException();
	}
	
	private void nextQueuedTask()
	{
		synchronized (queuedTasks)
		{
			for (Task task : queuedTasks.values())
			{
				CalculationServer server = ServerPool.getInstance().getFreeServer(task.taskType);
				if (server == null)
					throw new RuntimeException("Unbelievable! No free server for task type " + task.taskType);
				synchronized (server)
				{
					if (!server.busy)
					{
						server.busy = true;
						runningTasks.put(task.id, new CalculationThread(server, task));
						queuedTasks.remove(task.id);
						break;
					}
				}
			}
		}
		
	}
	
	public LocalTransport setEmbedded(boolean embedded) {
		this.embeddedTransport = embedded;
		return this;
	}
	
	public class CalculationThread extends Thread
	{
		CalculationServer server;
		Task task;
		
		public CalculationThread(CalculationServer server, Task task)
		{
			this.server = server;
			this.task = task;
			start();
		}
		
		public void run()
		{
			logger.info("Starting task " + task + " locally");
			server.calculateWrapper(task);
			logger.info("Finished task " + task + " locally");
			task.clearData();
			task.clearConfig();
			nextQueuedTask();
		}
	}
}
