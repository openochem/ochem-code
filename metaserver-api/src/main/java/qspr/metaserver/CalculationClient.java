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

import java.io.IOException;
import java.io.PrintWriter;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.CSTransport;
import qspr.metaserver.transport.Transport;
import qspr.metaserver.transport.TransportFactory;
import qspr.util.OverloadControl;
import qspr.util.StatusTracker;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.MemoryUtils;


@SuppressWarnings({"unchecked","rawtypes"})
public class CalculationClient
{
	public Transport transport = TransportFactory.create();
	public String sid = Long.valueOf(Math.round(Math.random()*1000000)).toString();
	public int interval = 1000;
	private String statusString = "idle";
	public Event statusChange = new Event(this);
	public Task currentTask;
	public PrintWriter out = new PrintWriter(System.out, true);
	private boolean tolerateIfMetaserverDown = false;

	public StatusTracker statusTracker;

	/**
	 * The ID of the user who is initiating the calculation.
	 */
	public String user;

	public Task calculateTask(Task task) throws Exception{
		return calculateTask(task, true);
	}

	public void setTolerateMetaserverDown(){
		tolerateIfMetaserverDown = true;
	}

	// Synchronous mode calculation
	public Task calculateTask(Task task, boolean alsoExecuteLocally) throws Exception
	{
		Task returnedTask;

		System.gc();
		int availableMemory = MemoryUtils.getCurrentMemoryFree();
		logger.info("Available memory in CalculationClient " + availableMemory + "MB required memory "+task.getMinRequiredMemory()+ "MB for "+task);

		// First trying to calculate task locally WITHOUT considering memory requirements
		// This is faster and less memory demanding than to post task, which involves serialization of objects
		if ((returnedTask = calculateTaskLocally(task)) != null)  
			return returnedTask;

		int taskId = postTask(task, alsoExecuteLocally);

		try
		{
			while (true)
			{
				returnedTask = getTask(taskId);
				if (returnedTask != null)
					return returnedTask;
				Thread.sleep(interval);
			}
		}
		catch (InterruptedException e)
		{
			// If we interrupt the thread - kill the task! We do care about performance / Midnighter
			logger.warn("[Client "+sid+"]: Killing task "+task+", because running thread has been interrupted");
			transport.executeCommand(new Command(Command.KILL_TASK, taskId).sid(sid));
			throw e;
		}
	}

	public Integer postTask(Task task) throws IOException, ClassNotFoundException, InterruptedException
	{
		return postTask(task, true);
	}

	public void setDeepSleepTime(int deepSleepTime)
	{
		if (transport instanceof CSTransport)
			((CSTransport) transport).deepSleepTime = deepSleepTime;
	}

	public List<Object[]> getTasksSummary() throws MalformedURLException, IOException, ClassNotFoundException {
		Command res = transport.executeCommand(new Command(Command.CL_GET_TASKS_SUMMARY, null).sid(sid));
		return (List<Object[]>) res.data;
	}

	//Alternate async concept - post a task, then return and pick it up
	public Integer postTask(Task task, boolean firstCalculateLocally) throws IOException, ClassNotFoundException, InterruptedException
	{
		task.setUser(user);

		long timeStart = Calendar.getInstance().getTimeInMillis();
		Command req = new Command(Command.CL_SUBMIT_TASK, task).sid(sid);
		Command res = null;
		boolean taskPosted = false;
		int attempts = 0;

		String md5 = task.getMD5();
		task.md5 = md5;
				
		if (md5 != null && task.isCachable())
		{
			firstCalculateLocally = false; // all cachable tasks are 
			Integer reference_id = getTaskByMD5(md5);
			if (reference_id != null)
			{
				task.id = reference_id;
				logger.info("[Client "+sid+"]: Task " + task + " has been found in task-level cache by md5("+md5+"), task size could have been " + task.getByteSize());
				out.flush();
				return reference_id;
			} else
				logger.info("[Client "+sid+"]: Task " + task + " has NOT been found in task-level cache by md5("+md5+"), proceeding with posting");
		}

		while (!taskPosted)
			try
		{
				if (transport instanceof CSTransport)
					res = ((CSTransport)transport).executeCommand(req, firstCalculateLocally);
				else
					res = transport.executeCommand(req);
				taskPosted = true;
		}
		catch (IOException e)
		{
			// Tolerate if the metaserver is down. But not more than 12 attempts, since this is a synchronious call! / Midnighter
			if (tolerateIfMetaserverDown && attempts++ < 12)
			{
				setStatus("Failed to post a task. Retrying in "+ attempts+" sec");
				Thread.sleep(1000*attempts);
			}
			else
				throw e;
		}
		Task postedTask = (Task)res.data;
		task.id = postedTask.id;
		logger.info("[Client "+sid+"]: Task " + task + " has been posted, it took "+(Calendar.getInstance().getTimeInMillis()-timeStart)/1000+" sec, task size was " + task.getByteSize());
		out.flush();
		return postedTask.id;
	}

	public Task getTask(Integer taskId) throws IOException, ClassNotFoundException
	{
		Command res;
		try
		{
			res = transport.executeCommand(new Command(Command.CL_QUERY_TASK, taskId).sid(sid));
		}
		catch (IOException e)
		{
			// 15 Feb 2010 Midnighter
			// It is now possible to tolerate the metaserver unavailability
			if (tolerateIfMetaserverDown)
			{
				logger.warn(Command.MEASERVER_DOWN);
				logger.warn("Error while connecting to the metaserver: " + e.getMessage());
				e.printStackTrace(out);
				return null;
			}
			else
				throw e;
		}

		if (res.id == Command.MS_UNKNOWN_TASK)
			throw new IOException("Metaserver lost our task "+taskId+"!");
		if (res != null)
		{
			Task mytask = currentTask = Task.fromCommand(res);
			if (mytask != null)
			{
				if (!Task.isAliveStatus(mytask.status))
				{
					logger.info("[Client "+sid+"]: Task "+mytask+" received");
					setStatus("Finished");
					return mytask;
				}

				if (statusString == null || !statusString.equals(mytask.getDetailedStatus()))
					setStatus(mytask.getDetailedStatus());
			}
		}

		return null;
	}

	public void setTaskPriority(int taskId, int priority) throws MalformedURLException, IOException, ClassNotFoundException
	{
		Task sendTask = new Task();
		sendTask.id = taskId;
		sendTask.setPriority(priority);

		try
		{
			transport.executeCommand(new Command(Command.CL_SET_PRIORITY, sendTask).sid(sid));
		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw new UserFriendlyException("Can't update the task priority while the metaserver is down. Please, retry later.");
		}
	}

	public Task getReadyTask(String taskType) throws Exception
	{
		Command res = transport.executeCommand(new Command(Command.CL_QUERY_READY_TASK, taskType).sid(sid));
		return currentTask = Task.fromCommand(res);
	}

	/**
	 * Check if the metaserver is already calculating a task with given MD5 hash.
	 * Used for avoiding of posting duplicate tasks
	 * @throws IOException 
	 */
	private Integer getTaskByMD5(String md5) throws IOException, ClassNotFoundException, InterruptedException //Should it be blocking, like it is now? Only used in blocking "post task", so I guess it won't harm
	{
		OverloadControl oc = new OverloadControl("Metaserver", 1, 10000);
		Command req = new Command(Command.CL_GET_TASK_BY_MD5, md5).sid(sid);
		Command res = null;
		boolean gotResult = false;

		while (!gotResult)
		{
			try
			{
				res = transport.executeCommand(req);
				if (res == null)
					return null;
				gotResult = true;
			} catch (IOException e)
			{
				if (!tolerateIfMetaserverDown || oc.getNumMaxTimeout() > 12)
					throw e;
				setStatus("Failed to query task by md5. Retrying in " + (oc.getCurrentTimeout() / 1000) + " sec");
				oc.relax(e);
			}
		}
		return (Integer) res.data;
	}


	public Task getTaskStatus(Integer taskId, boolean retrieveResult) throws IOException, ClassNotFoundException
	{

		int cmdCode = retrieveResult ? Command.CL_QUERY_TASK : Command.CL_QUERY_TASK_STATUS; 
		Command res;
		try
		{
			res = transport.executeCommand(new Command(cmdCode, taskId).sid(sid));
		}
		catch (IOException e)
		{
			if (tolerateIfMetaserverDown)
			{
				logger.warn(Command.MEASERVER_DOWN);
				logger.warn(e.getMessage());

				// Return a dummy task
				Task task = new Task();
				task.status = Task.INIT;
				task.setDetailedStatus(Command.MEASERVER_DOWN);
				return task;
			}
			else
				throw e;
		}
		if (res.id == Command.MS_UNKNOWN_TASK)
			throw new IOException("Metaserver lost our task "+taskId+"!");

		Task mytask = Task.fromCommand(res);
		return mytask;
	}

	public List<Task> getTaskStatuses(List<Integer> taskIds) throws Exception
	{
		// A method to request statuses of multiple tasks with a single query to metaserver
		// Not supported for local tasks!

		ArrayList<Task> tasks = new ArrayList<Task>();
		if(taskIds.isEmpty()) return tasks;
		Command res;

		while (taskIds.size() > 10)
		{
			tasks.addAll(getTaskStatuses(taskIds.subList(0, 10)));
			taskIds = taskIds.subList(10, taskIds.size());
		}

		try
		{
			long time = Calendar.getInstance().getTimeInMillis();

			if (!(taskIds instanceof ArrayList))
			{
				List<Integer> newList = new ArrayList<Integer>();
				newList.addAll(taskIds);
				taskIds = newList;
			}

			res = transport.executeCommand(new Command(Command.CL_QUERY_TASK_STATUS, (ArrayList<Integer>)taskIds).sid(sid));
			logger.info("Fetched statuses of " + taskIds.size() + " tasks in " + (Calendar.getInstance().getTimeInMillis() - time) + " ms.");
		}
		catch (IOException e)
		{
			if (tolerateIfMetaserverDown)
			{
				logger.warn(Command.MEASERVER_DOWN);
				logger.warn(e.getMessage());
				for (int i = 0; i < taskIds.size(); i++)
				{
					// Return a dummy task
					Task task = new Task();
					task.id = taskIds.get(i);
					task.status = Task.INIT;
					task.setDetailedStatus(Command.MEASERVER_DOWN);
					tasks.add(task);
				}
			}
			else
				throw e;

			return tasks;
		}

		if (tasks.isEmpty())
			return (List<Task>)res.data;
		else
		{
			tasks.addAll((List<Task>)res.data);
			return tasks;
		}
	}

	/**
	 * Retrieve the list of the task types supported by at least one available calculation server.
	 * @return
	 */
	public Set<String> getSupportedTaskTypes() throws MalformedURLException, IOException, ClassNotFoundException
	{
		Command res = transport.executeCommand(new Command(Command.CL_GET_SUPPORTED_TASKS, null).sid(sid));
		return (Set<String>) res.data;
	}

	public void killTask(Integer taskId) throws IOException, ClassNotFoundException
	{
		logger.info("Sending kill for "+taskId);
		try
		{
			transport.executeCommand(new Command(Command.KILL_TASK, taskId).sid(sid));
		}
		catch (IOException e)
		{
			if (tolerateIfMetaserverDown)
				logger.warn("Failed to kill a task " + taskId);
			else
				throw e;
		}
	}

	public void deleteTask(Integer taskId) throws Exception
	{
		logger.info("Sending delete for "+taskId);
		transport.executeCommand(new Command(Command.CL_DELETE_TASK, taskId).sid(sid));
	}

	public Integer getParentId(Integer taskId) throws MalformedURLException, IOException, ClassNotFoundException, InterruptedException
	{
		OverloadControl oc = new OverloadControl("Metaserver", 1, 10000);
		Command req = new Command(Command.CL_GET_PARENT_ID, taskId).sid(sid);
		Command res = null;
		boolean gotResult = false;

		while (!gotResult)
		{
			try
			{
				res = transport.executeCommand(req);
				if (res == null)
					return null;
				gotResult = true;
			} catch (IOException e)
			{
				if (oc.getNumMaxTimeout() > 12)
					throw e;
				setStatus("Failed to query task by md5. Retrying in " + (oc.getCurrentTimeout() / 1000) + " sec");
				oc.relax(e);
			}
		}
		return (Integer) res.data;
	}

	public void setTaskParent(Integer childID, Integer parentID) throws MalformedURLException, IOException, ClassNotFoundException
	{
		transport.executeCommand(new Command(Command.CL_SET_PARENT, new Integer[]{parentID, childID}).sid(sid));	
	}

	public String getStatus()
	{
		return statusString;
	}

	public void setStatus(String statusString)
	{
		if (statusString == null || !statusString.equals(this.statusString))
		{
			this.statusString = statusString;
			statusChange.fire();
			if (statusTracker != null)
				statusTracker.set(statusString);
		}
	}

	public CalculationClient setSid(String sid)
	{
		this.sid = sid;
		return this;
	}

	private Task calculateTaskLocally(Task task) throws IOException
	{
		// There is a local pool of servers. Use a local server if possible / Midnighter
		final CalculationServer server = 
				ServerPool.getInstance().getFreeServer(task.taskType);

		if (server != null)
		{
			server.statusChange.addListener(new EventListener<String>()
					{
				public void onEvent(Event event, String arg) 
				{
					CalculationClient.this.setStatus(server.getStatus());
				}
					});
			task.id = -(int)Math.round(Math.random()*1000);
			task.calcServerId = ServerPool.getInstance().sid;
			server.calculateWrapper(task);
			return task;
		}

		if(server == null && task.isLocalTask())throw new IOException("The task marked as local cannot be executed locally "+task);

		return null;
	}

	public CalculationClient()
	{
	}

	public CalculationClient setStatusTracker(StatusTracker status) {
		this.statusTracker = status;
		return this;
	}

	public CalculationClient(String sid)
	{
		this.sid = sid;
	}

	public CalculationClient(String sid, String user)
	{
		this.sid = sid;
		this.user = user;
	}

	private static Logger logger = LogManager.getLogger(CalculationClient.class);

}

