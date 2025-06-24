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

import java.io.ByteArrayOutputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.time.LocalTime;
import java.time.temporal.ChronoUnit;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.flywaydb.core.Flyway;
import org.hibernate.Transaction;
import org.hibernate.cfg.Configuration;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.type.IntegerType;
import org.hibernate.type.StringType;

import qspr.metaserver.frontend.SimpleController;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.ServerInfo;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.DataReferenceCleaner;
import qspr.metaserver.transport.NoSQLCleanerThread;
import qspr.metaserver.transport.NoSqlTransport;
import qspr.updateserver.UpdateServer;
import qspr.util.DaemonThread;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.config.ConfigurableClass;
import com.eadmet.utils.config.ConfigurableProperty;
import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;

@ConfigurableClass(name = "metaserver", comment = "Meta server configs")
@SuppressWarnings("unchecked")
public class MetaServer extends AbstractServer implements DataReferenceCleaner
{
	private static transient final Logger logger = LogManager.getLogger(MetaServer.class);
	public static final String METASERVER = "metaserver";
	final static double PRIORITY_DECREMENT_FOR_RESTART = 0.1;
	public static final String IMMEDIATE_FAILURE = "FAILURE:";

	// Couple of non-automatic cleaning commands
	// delete Task from Task left join ochem.PendingTask using(task_id) where Task.status="ready" and PendingTask.task_id is NULL and parent_task_id is null;
	// delete from Task where status="error" or status="kill";
	// delete from Task where  parent_task_id is not null;
	// delete from Task where task_type  not in ("Workflow","Consensus"); -- only if there are no other running tasks
	// update ochem.PendingTask set published=0 where article_id is  NULL;
	// select count(*) from ochem.PendingTask left join metaserver.Task using(task_id) where Task.task_id is NULL and published =0 ;

	// Update of published Baskets to make them publicly available:
	// update Basket join Model  on (training_set_id=basket_id)  set user_id=1 where published=1 and user_id!=1 and (published_id<630 or published_id>635);
	// update Basket join Model  on (validation_set_id=basket_id)  set user_id=1 where published=1 and user_id!=1 ;
	// update Basket join ValidationSet using (basket_id) join Model using (model_id) set user_id = 1  where published=1 and user_id!=1;

	@ConfigurableProperty(name = "maindb")
	public static Map<String, String> mainDbOverrides;

	/**
	 * Defines whether this is a mirror installation of the Metaserver (e.g., developer machine)
	 */
	@ConfigurableProperty(name = "mirror")
	public static boolean mirror;

	public static final int METASERVER_TASK_OVERDUE = 60 * 60 * 1000;
	public static final int METASERVER_TASK_NOT_RESPONING = 5 * 60 * 1000; // Time to consider the server as non-responding (gray it out)
	public static final int RESUBMIT_FAILED_TASK_ATTEMPTS = 5;
	public static final int LIMITCHILDREN = 200; // to prevent very long queries

	/**
	 * How long can exclusive tasks wait for "their" server?
	 */
	private static final int EXCLUSIVE_TASKS_TIMEOUT = 200 * 1000;
	private static final String SEPARATOR = " -- ";
	private static final String STOP_REQUSTED = " manual stop requested ";

	public long startupTime;

	/**
	 * The map of queued tasks.
	 * Only one task per task type is stored in the queue.
	 */
	public TaskQueue taskQueue;


	/**
	 * The map of currently assigned tasks.
	 * Maps the server id (SID) to the task being calculated by this server
	 */
	public Map<String, Task> assignedTasks; // "MLRA_AT_Harthof" -> Task()

	/**
	 * The map of online peers (both servers and clients).
	 * Maps the sender ID (SID) to the OnlinePeer object
	 */
	public Map<String, OnlinePeer> onlinePeers;

	/**
	 * A map of timestamps when a task type appeared on the metaserver.
	 * Used to detect "dead" task types that do not have any supporting servers.
	 */
	public Map<String, Long> taskTypePing = new HashMap<String, Long>();


	public String currentVersion;
	public boolean terminatorMode;
	public String currentCommand;
	long lockStartedTime;

	/**
	 * A flag to indicate that the metaserver should stop assigning new tasks.
	 * This is useful to control overloads.
	 */
	public boolean stopTaskAssignment = false;

	/**
	 * This is used to forcely send all task of a particular type to a particular server.
	 * Useful for "global" debugging of a task type.
	 */
	public Map<String, String> forcedTaskBindings = new HashMap<String, String>();

	private static MetaServer instance;

	protected AtomicInteger newTasksCounter = new AtomicInteger(0);
	protected AtomicInteger completedTasksCounter = new AtomicInteger(0);
	protected AtomicInteger errorsCounter = new AtomicInteger(0);

	@Override
	protected Configuration configureDbConnection(Configuration hibernateConf) throws Exception
	{
		if (mainDbOverrides != null)
			for (String key : mainDbOverrides.keySet())
				hibernateConf.setProperty(key, mainDbOverrides.get(key));
		return hibernateConf;
	}

	public Command executeCommand(Command request)
	{
		Command response;
		long time = Calendar.getInstance().getTimeInMillis();
		long actualTime;
		debug("         -th queue number execution <" + request + "> by " + request.senderId + " " + " with id " + Thread.currentThread().getId());
		synchronized (this)
		{
			lockStartedTime = Calendar.getInstance().getTimeInMillis();
			currentCommand = request + "(" + request.senderId + ")";
			debug("         0th queue number execution <" + request + "> by " + request.senderId + " " + " with id " + Thread.currentThread().getId());
			Transaction tx = session().beginTransaction();

			try
			{
				synchronized (this)
				{
					debug("         1th queue number execution <" + request + "> by " + request.senderId + " " + " with id " + Thread.currentThread().getId());
					actualTime = Calendar.getInstance().getTimeInMillis();
					response = executeCommandInternal(request);
					session().flush();
					tx.commit();
					actualTime = Calendar.getInstance().getTimeInMillis() - actualTime;
					debug("         2th queue number execution <" + request + "> by " + request.senderId + " " + " with id " + Thread.currentThread().getId());
				}
			} catch (Exception e)
			{
				tx.rollback();
				throw new RuntimeException(e);
			}
			debug("         3th queue number execution <" + request + "> by " + request.senderId + " " + " with id " + Thread.currentThread().getId());
			currentCommand = null;
		}

		/*Try to determine packets size*/
		int reqSize = 0;
		int resSize = 0;
		try
		{
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			ObjectOutputStream oos = new ObjectOutputStream(baos);
			oos.writeObject(request);
			oos.flush();
			reqSize = baos.size();
			oos.reset();
			baos.reset();
			oos.writeObject(response);
			oos.flush();
			resSize = baos.size();
			oos.close();
		} catch (Exception e)
		{
			e.printStackTrace();
		}
		if (request.id != Command.CS_REGISTER || resSize > 100) // Do not show empty server registrations not to trush the logs
			logger.info(request + "(" + request.senderId + "):" + reqSize + " - " + response + ":" + resSize + "  "
					+ (Calendar.getInstance().getTimeInMillis() - time) + "ms. (" + actualTime + " ms.)");
		return response;
	}

	protected Command executeCommandInternal(Command command)
	{
		// Do not accept anonymous commands
		if (command.senderId == null)
			return new Command(Command.DENIAL, null);

		OnlinePeer peer = addPeer(command.senderId);

		try
		{
			switch (command.id.intValue())
			{
			case Command.CL_SUBMIT_TASK:
				return clSubmitTask(command, peer);

			case Command.CL_SET_PRIORITY:
				return clSetPriority(command, peer);

			case Command.CL_GET_SUPPORTED_TASKS:
				return clGetSupportedTasks(command, peer);

			case Command.CL_GET_FAILED_TASKS:
				return clGetFailedTasks(command, peer);

			case Command.CL_SET_PARENT:
				return clSetParent(command, peer);

			case Command.CL_REGISTER_ADMIN_IP:
				return clRegisterAdminIP(command, peer);

			case Command.CL_GET_TASK_BY_MD5:
				return clGetTaskByMD5(command, peer);

			case Command.CL_GET_PARENT_ID:
				List<Number[]> rows =  session().createSQLQuery("select parent_task_id from Task where task_id=:id").setParameter("id", command.getDataInteger()).list();
				return rows == null || rows.size() == 0 ? null : new Command(Command.MS_TASKFOUND, rows.get(0));

			case Command.CL_QUERY_TASK:
				return clQueryTask(command, peer);

			case Command.CL_QUERY_TASK_STATUS:
				return clQueryTaskStatus(command, peer);

			case Command.CL_QUERY_READY_TASK:
				return clQueryReadyTask(command, peer);

			case Command.CL_GET_TASKS_SUMMARY:
				return clGetTasksSummary(command, peer);

			case Command.CL_DELETE_TASK:
				logger.info("Deleting task " + command.data + " as requested by " + command.senderId);
				session().createQuery("delete from Task where id=:id").setInteger("id", command.getDataInteger()).executeUpdate();
				removeTaskFromCache(command.getDataInteger());
				deleteChildrenTasks(command.getDataInteger(), false, false);
				return null;
			case Command.CL_DELETE_CHILDREN:
				logger.info("Deleting children of " + command.data + " as requested by " + command.senderId);
				deleteChildrenTasks(command.getDataInteger(), false, false);
				return null;
			case Command.CS_REGISTER:
				return registerServer(command);

			case Command.CS_TASK_CALCULATED:
				return csTaskCalculated(command, peer);

			case Command.CS_QUERYINFO:
				return csQueryInfo(command, peer);

			case Command.KILL_TASK:
				return killTask(command, peer);

			default:
				return null;
			}
		}
		finally
		{
			registerServerConnection(peer);
		}
	}

	private Command clGetTaskByMD5(Command command, OnlinePeer peer)
	{
		String md5 = (String) command.data;
		Integer taskId = getTaskByMD5(md5);

		if (taskId == null)
		{
			logger.info("NOT found a task by MD5 " + md5);
			return null;
		}
		else
		{
			logger.info("Found a task by MD5 + " + md5 + ". Task ID is " + taskId);
			return new Command(Command.MS_TASKFOUND, taskId);
		}
	}

	private Command clSetParent(Command command, OnlinePeer peer)
	{
		Integer[] taskIDs = (Integer[]) command.data;
		session().createSQLQuery("update Task set parent_task_id=:parent where task_id=:child")
		.setInteger("parent", taskIDs[0])
		.setInteger("child", taskIDs[1]).executeUpdate();

		Task task = getTaskFromCache(taskIDs[1]);
		if (task != null)
			task.parentTaskId = taskIDs[0];

		return null;
	}

	private Integer getTaskByMD5(String md5)
	{
		List<Object[]> rows =  session().createSQLQuery("select status, task_id from Task where task_md5=:md5  and status != \\\"error\\\" and  status != \\\"kill\\\"").setParameter("md5", md5).list();
		Integer taskId = null;
		for (Object[] row : rows)
		{
			String status = (String) row[0];
			Integer id = ((Number) row[1]).intValue();
			if (Task.READY.equals(status))
				taskId = id;
			else if (Task.isAliveStatus(status))
				if (taskId == null)
					taskId = id;
		}

		if (taskId != null)
		{
			// We got a task!
			Task task = getTaskFromCache(taskId);
			if (task != null)
				task.referenceCount++;
			session().createSQLQuery("update Task set ref_count=ref_count+1 where task_id=:id").setParameter("id", taskId).executeUpdate();
		}

		return taskId;
	} 

	boolean isToBeResubmitted(Task task){
		return isParentAlive(task.parentTaskId) && // parent is alive or it is really small task
				(task.resubmitted == null || task.resubmitted < RESUBMIT_FAILED_TASK_ATTEMPTS) && // resubmitted less than allowed number
				!Task.KILL.equals(task.status); // was not killed
	}

	boolean isParentAlive(Integer parentTaskId){
		if(parentTaskId == null) return false; // there is no parent!
		String parentStatus = (String) session().createSQLQuery("select status from Task where task_id=:taskId")
				.setParameter("taskId", parentTaskId).uniqueResult();
		return Task.isAliveStatus(parentStatus);
	}

	/**
	 * Process the new task submission
	 * @return response command
	 */
	private Command clSubmitTask(Command command, OnlinePeer peer)
	{
		Task task = (Task) command.data;

		task.referenceCount = 1;
		task.lastAccess = new Timestamp(Calendar.getInstance().getTimeInMillis());
		task.status = Task.INIT;

		double initialPriority = task.priority;

		if(task.parentTaskId != null) {
			Task parent = (Task) session().get(Task.class, task.parentTaskId);
			initialPriority = parent != null ? parent.priority : initialPriority;
		}

		task.correctTaskPriority(initialPriority);
		task.time = new Timestamp(Calendar.getInstance().getTimeInMillis());

		if (task.md5 != null && task.isCachable())
		{
			Integer taskId = getTaskByMD5(task.md5);
			if (taskId != null)
			{
				Task resTask = new Task();
				resTask.id = taskId;
				resTask.status = task.status;
				resTask.time = task.time;
				return new Command(Command.MS_TASK_STATUS, resTask);
			}
		}

		if (task.parentTaskId !=null && !isParentAlive(task.parentTaskId))
		{
			logger.info("Task " + task + " is a child of a killed task and, therefore, it has been disregarded");
			task.status = Task.KILL;
			task.setDetailedStatus("Task disregarded as a child of an unknown/killed task");
		}

		task.client = command.senderId;
		if (task.taskName != null)
			task.taskName = task.taskName.substring(0, Math.min(task.taskName.length(), 99));


		checkTaskSize(task);

		// Update user from the parent task
		if (task.getUser() == null && task.parentTaskId != null)
			task.setUser(getTaskUser(task.parentTaskId));

		session().save(task);

		// Put this task into the cache, if there is no task or if there is a task with a lower priority
		taskQueue.putTask(task);

		peer.status = "submitted a task";
		peer.isClient = true;
		peer.currentTask = task;

		// Dont send task data back
		Task resTask = new Task();
		resTask.id = task.id;
		resTask.status = task.status;
		resTask.time = task.time;

		newTasksCounter.addAndGet(1);

		return new Command(Command.MS_TASK_STATUS, resTask);
	}

	private Command clGetTasksSummary(Command command, OnlinePeer peer) {
		List<Object[]> rows = session().createSQLQuery("select task_type, status, count(*) cnt from Task group by task_type, status")
				.addScalar("task_type", StringType.INSTANCE)
				.addScalar("status", StringType.INSTANCE)
				.addScalar("cnt", IntegerType.INSTANCE)
				.list();

		return new Command(Command.MS_OK, (Serializable) rows);
	}

	/**
	 * Add an IP address to the list of authorized addresses
	 */
	private Command clRegisterAdminIP(Command command, OnlinePeer peer)
	{
		String newAddress = (String) command.data;
		Calendar newExpiry = Calendar.getInstance();
		newExpiry.add(Calendar.HOUR, 1);

		logger.info("Current admin IPs:");
		SimpleDateFormat simdf = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
		boolean updatedExisting = false;

		for (Object[] tuple : allowedAddresses)
		{
			String address = (String) tuple[0];
			Date expiry = (Date) tuple[1];
			if (newAddress.startsWith(address))
			{
				updatedExisting = true;
				if (expiry.before(newExpiry.getTime()))
				{
					tuple[1] = newExpiry.getTime();
					logger.info(address + " expires " + simdf.format(newExpiry.getTime()) + " (new address matched, expiry updated)");
				}
				else
					logger.info(address + " expires " + simdf.format(expiry) + " (new address matched)");
			}
			else
				logger.info(address + " expires " + simdf.format(expiry));
		}

		if (!updatedExisting)
		{
			Object[] tuple = { newAddress, newExpiry.getTime() };
			allowedAddresses.add(tuple);
			logger.info("Registered " + newAddress + " as admin IP, expires " + simdf.format(newExpiry.getTime()));
		}

		return null;
	}

	private Command clQueryTaskStatus(Command command, OnlinePeer peer)
	{
		List<Number> requestedTaskIds  = new ArrayList<Number>();
		ArrayList<Task> tasks = new ArrayList<Task>();

		if (command.data instanceof List)
			requestedTaskIds = (List<Number>) command.data;
		else
			requestedTaskIds.add((Number)command.data);

		if (requestedTaskIds == null || requestedTaskIds.isEmpty())
			return new Command(Command.MS_TASK_STATUS, tasks);

		List<Object[]> rows = session().createSQLQuery("select task_id, status, detailed_status, calc_server_id from Task where task_id in (:ids)")
				.setParameterList("ids", requestedTaskIds).list();

		// Update the last_access time
		session().createSQLQuery("update Task set last_access=:now where task_id in (:ids)")
		.setParameterList("ids", requestedTaskIds)
		.setParameter("now", Calendar.getInstance().getTime())
		.executeUpdate();

		for (int i = 0; i < requestedTaskIds.size(); i++)
		{
			Task task = new Task();
			task.id = requestedTaskIds.get(i).intValue();
			task.setError(Command.UNKNOWN_TASK + task.id);

			for (Object[] row : rows)
				if (task.id == ((Number)row[0]).intValue())
				{
					task.status = (String) row[1];
					task.setDetailedStatus((String) row[2]);
					task.calcServerId = (String) row[3];
					if (Task.INIT.equals(task.status))
						task.setDetailedStatus("Waiting for a free server");
				}

			tasks.add(task);
		}

		if(command.data instanceof List)
			return new Command(Command.MS_TASK_STATUS, tasks);
		else
			return new Command(Command.MS_TASK_STATUS, tasks.get(0));
	}

	protected Command csTaskCalculated(Command command, OnlinePeer peer)
	{
		Task task = (Task) command.data;

		Task storedTask = (Task) session().get(Task.class, task.id);
		if (storedTask == null)
		{
			logger.info(String.format("%s reports an unknown or deleted task %s ", command.senderId, task));
			return null;
		}

		if (!"ready".equals(task.status) && !command.senderId.equalsIgnoreCase(storedTask.calcServerId))
		{
			// Wait a moment.. This is not your task. Ignore. (Possible, when task was resubmitted to somebody else)
			logger.info(String.format("%s reported the task %s officially calculated by %s", command.senderId, task, storedTask.calcServerId));
			return null;
		}

		if ("ready".equals(storedTask.status))
			return null;

		deleteChildrenTasks(storedTask.id, true, false);

		storedTask.timeCompleted = new Timestamp(Calendar.getInstance().getTimeInMillis());
		storedTask.setDetailedStatus(task.getDetailedStatus());
		storedTask.updateResult(task);
		storedTask.status = task.status;
		storedTask.peakMemoryUsage = task.peakMemoryUsage;
		checkTaskSize(storedTask);
		peer.currentTask = null;

		completedTasksCounter.addAndGet(1);
		assignedTasks.remove(storedTask.calcServerId);

		// If this is a part of another task and the failure is temporary, resubmit the task
		if (Task.ERROR.equals(storedTask.status))
		{
			logger.info("Task " + task + " has finished with an error: " + task.getDetailedStatus());
			errorsCounter.addAndGet(1);
			if (!task.getDetailedStatus().startsWith(IMMEDIATE_FAILURE) && isToBeResubmitted(storedTask))
			{	
				resubmitTask(storedTask);
				return null;
			}
		}

		logger.info("Task " + storedTask + " is finished.");

		if (!Task.Workflow.equals(storedTask.taskType)) // we zero data for all child tasks to save memory when sending them to the server -- added by IVT 17.08.2011
			storedTask.clearData();
		peer.status = "has finished task " + storedTask;
		session().saveOrUpdate(storedTask);
		session().save(new ArchivedTask(storedTask)); // Store the task in the archive table

		return null;
	}

	protected Command csQueryInfo(Command command, OnlinePeer peer)
	{
		Task task;
		ServerInfo receivedServerInfo = (ServerInfo) command.data;
		if (receivedServerInfo.configurationXml == null && peer.serverInfo != null)
			receivedServerInfo.configurationXml = peer.serverInfo.configurationXml;
		ServerInfo serverInfo = peer.serverInfo = receivedServerInfo;

		synchronized (assignedTasks)
		{
			if ((task = assignedTasks.get(command.senderId)) == null)
			{
				logger.info("Server " + command.senderId + "(" + peer.serverInfo.random + ") is calculating an unknown or deleted task, the status was "
						+ peer.serverInfo.status);
				peer.status = "Received a kill task signal";
				peer.currentTask = null;
				return new Command(Command.MS_UNKNOWN_TASK, null);
			}
		}

		if (Task.KILL.equals(task.status))
		{
			task.status = Task.KILLED;
			session().createQuery("update Task set status = :status where id = :id").setString("status", task.status).setInteger("id", task.id)
			.executeUpdate();
			peer.status = "Received a kill task signal";
			peer.currentTask = null;
			return new Command(Command.KILL_TASK, task.id);
		}

		peer.currentTask = task;

		if (Task.STOP.equals(task.status))
		{
			updateDetailedStatus(peer, (String) serverInfo.status);
			return new Command(Command.MS_STOP, task.id);
		}

		if (serverInfo.status != null)
		{
			String newStatus = (String) serverInfo.status;
			if(newStatus.length() > Task.MAX_STATUS_LENGTH) newStatus = newStatus.substring(0, Task.MAX_STATUS_LENGTH);
			if(!task.getDetailedStatus().startsWith(newStatus)) // something new!
				updateDetailedStatus(peer, newStatus);
		}
		else
			peer.status = null;

		if (peer.logsRequested)
		{
			peer.logsRequested = false;
			return new Command(Command.MS_GET_LOGS, null);
		}

		return null;
	}

	void updateDetailedStatus(OnlinePeer peer, String detailedStatus) {
		if(detailedStatus.length() > Task.MAX_STATUS_LENGTH) detailedStatus = detailedStatus.substring(0, Task.MAX_STATUS_LENGTH);
		detailedStatus += SEPARATOR + (Task.STOP.equals(peer.currentTask.status)?STOP_REQUSTED:"") + LocalTime.now().truncatedTo(ChronoUnit.MINUTES).toString();
		peer.currentTask.setDetailedStatus(detailedStatus);
		peer.status = detailedStatus;
		session().createQuery("update Task set detailedStatus = :dstatus where id = :id").setString("dstatus", detailedStatus)
		.setInteger("id", peer.currentTask.id).executeUpdate();
	}

	protected Command clQueryTask(Command command, OnlinePeer peer)
	{
		Integer taskId = command.getDataInteger();

		Task task = (Task) session().get(Task.class, taskId);

		if (task == null)
			return new Command(Command.MS_UNKNOWN_TASK, null);
		task.lastAccess = new Timestamp(Calendar.getInstance().getTimeInMillis());
		session().save(task);

		peer.currentTask = task;
		peer.isClient = true;

		if (Task.INIT.equals(task.status))
			task.setDetailedStatus("Waiting for a free server");
		if (isTaskOverdue(task))
			onTaskOverdue(task);

		if (task.isReady())
			peer.status = "received task results";
		else
			peer.status = "waiting";

		session().evict(task);
		// Dont send task data back
		task.clearData();

		return new Command(Command.MS_TASK_STATUS, task);
	}

	/**
	 * Special command for fragment-calculating tasks.
	 * Main idea is to fetch any ready task of specific type
	 * and remove it from MetaServer Tasks table completely
	 * NoS & Rob 09.02.10

			// Due to problem with dissappearing tasks now
			// we're not deleting them, but changing status
			// NoS 15.03.10
	 * @return
	 */
	protected Command clQueryReadyTask(Command command, OnlinePeer peer)
	{

		peer.isClient = true;
		String taskType = (String) command.data;

		Task task;
		Task readyTask = getTask("ready", taskType);
		if (readyTask == null)
			return new Command(Command.MS_UNKNOWN_TASK, null);
		else
			task = readyTask;

		peer.currentTask = task;

		task.status = Task.FETCHED;
		session().saveOrUpdate(task);
		session().flush();
		session().evict(task);

		// Dont send task data back
		task.clearData();

		return new Command(Command.MS_TASK_STATUS, task);
	}

	/**
	 * Get the list of task types non-supported by at least one server
	 */
	protected Command clGetFailedTasks(Command command, OnlinePeer peer)
	{
		HashSet<String> failedTasks = new HashSet<String>();
		for (OnlinePeer onlinePeer : onlinePeers.values())
		{

			if (onlinePeer.isClient || onlinePeer.getPing() > MetaServer.METASERVER_TASK_NOT_RESPONING)
				continue;
			if (onlinePeer.serverInfo != null && onlinePeer.serverInfo.failures != null)
				failedTasks.addAll(onlinePeer.serverInfo.failures);
		}
		return new Command(Command.MS_OK, failedTasks);
	}

	/**
	 * Get the list of task types supported by at least one server
	 */
	protected Command clGetSupportedTasks(Command command, OnlinePeer peer)
	{
		HashSet<String> supportedTasks = new HashSet<String>();
		for (OnlinePeer onlinePeer : onlinePeers.values())
		{

			if (onlinePeer.isClient || onlinePeer.getPing() > MetaServer.METASERVER_TASK_NOT_RESPONING)
				continue;
			if (onlinePeer.serverInfo != null && onlinePeer.serverInfo.supportedTaskTypes != null)
				supportedTasks.addAll(onlinePeer.serverInfo.supportedTaskTypes);
		}
		return new Command(Command.MS_OK, supportedTasks);
	}

	private Command killTask(Command command, OnlinePeer peer)
	{
		Integer taskId = command.getDataInteger();
		Task task = (Task) session().get(Task.class, taskId);
		logger.info("" + task + " killed by " + command.senderId);
		if (task == null)
			return new Command(Command.MS_UNKNOWN_TASK, null);
		if (Task.isAliveStatus(task.status))
		{
			if (!task.hasMultipleReferences())
			{
				task.status = Task.KILL;
				peer.status = "requested a murder of a task";
				task.setDetailedStatus("Killed by the client");
			}
			else
				task.scheduledKill = true;
		}
		session().saveOrUpdate(task);
		deleteChildrenTasks(task.id, true, false);
		if (removeTaskFromCache(taskId))
			queueNextTask(task);
		return null;
	}


	private Command clSetPriority(Command command, OnlinePeer peer)
	{
		Task task = (Task) command.data;
		task.correctTaskPriority(task.priority);
		session().createSQLQuery("update Task set priority=:priority where (task_id=:taskId) or (parent_task_id=:taskId)") // Also need to consider children of children.. Requires more SQL requests.
		.setParameter("priority", task.priority).setParameter("taskId", task.id).executeUpdate();
		reloadCache(false, true); // Alternatively, we could use a more optimal but a more complex way to update the cache
		return null;	
	}

	/**
	 * Register the event of the server connection
	 */
	private void registerServerConnection(OnlinePeer peer)
	{
		if (peer.serverInfo != null && peer.serverInfo.supportedTaskTypes != null)
		{
			long time = Calendar.getInstance().getTimeInMillis();
			for (String taskType : peer.serverInfo.supportedTaskTypes)
				taskTypePing.put(taskType, time);
		}
	}

	private String getTaskUser(int taskId)
	{
		Task cachedTask = getTaskFromCache(taskId);
		if (cachedTask != null)
			return cachedTask.getUser();
		return (String) session().createSQLQuery("select user from Task where task_id=:id").setParameter("id", taskId).uniqueResult();
	}

	/**
	 * Change task status
	 * @param task
	 */

	private void changeTaskStatus(Task task){
		killChildrenTasks(task.id);
		removeTaskFromCache(task.id);
		if(task.calcServerId!=null)assignedTasks.remove(task.calcServerId);
		session().saveOrUpdate(task); // this will save new status
		reloadCache(true, true);
	}

	/**
	 *  Tasks will still remain at the servr
	 * @param parentId
	 */
	private void killChildrenTasks(Integer parentId)
	{
		if (parentId == null)
			return;

		List<Number> taskIds = session().createCriteria(Task.class).add(Restrictions.eq("parentTaskId", parentId)).setProjection(Projections.id()).list();

		if (taskIds.size() > 0)
		{
			session()
			.createQuery(
					"update Task set result=null, data=null, status='kill', detailedStatus='Killed by metaserver as a child of a killed task' where (ref_count <= 1) and id in (:id) and (status = 'init' or status='assigned')")
			.setParameterList("id", taskIds).executeUpdate();

			session()
			.createQuery(
					"update Task set scheduled_kill=1 where (ref_count > 1) and id in (:id)")
			.setParameterList("id", taskIds).executeUpdate();

			for (Number id : taskIds)
			{
				logger.info("Killing " + id + "." + parentId + " as a child of " + parentId);
				removeTaskFromCache(id.intValue());
				killChildrenTasks(id.intValue());
			}

			session().flush();
		}

	}

	private void deleteChildrenTasks(Integer parentId, boolean keepErrors, boolean cleanUpOnly) {}
	/*
	private void deleteChildrenTasks(Integer parentId, boolean keepErrors, boolean cleanUpOnly)
	{
		if (parentId == null)
			return;
		List<Number> taskIds = session().createCriteria(Task.class).add(Restrictions.eq("parentTaskId", parentId)).setProjection(Projections.id()).list();
		for (Number id : taskIds)
		{
			logger.info("Deleting " + id + " as a child of " + parentId);
			if (cleanUpOnly)
			{
				// Cleanup a task
				session().createQuery("update Task set result=null, data=null where ref_count <= 1 and id=:id" + (keepErrors ? "  and status<>'error'" : "")).setInteger("id", id.intValue())
				.executeUpdate();
				// Kill the task, if still active
				session().createQuery("update Task set status='kill' where ref_count <= 1 and id=:id and (taskType='init' or taskType='assigned')").setInteger("id", id.intValue())
				.executeUpdate();
			}
			else
				session().createQuery("delete from Task where ref_count <= 1 and id=:id" + (keepErrors ? "  and status<>'error'" : "")).setInteger("id", id.intValue()).executeUpdate();

			session().createQuery("update Task set scheduled_kill=1 where id=:id and ref_count > 1").setInteger("id", id.intValue())
			.executeUpdate();

			removeTaskFromCache(id.intValue());
			deleteChildrenTasks(id.intValue(), keepErrors, cleanUpOnly);
		}
		session().flush();
	}
	 */
	private boolean removeTaskFromCache(int taskId)
	{
		String foundKey = null;

		for (String key : assignedTasks.keySet())
			if (assignedTasks.get(key).id.equals(taskId))
				foundKey = key;
		if (foundKey != null)
			assignedTasks.remove(foundKey);

		return taskQueue.removeTaskById(taskId);
	}

	// Find the task with the given ID in our cache
	private Task getTaskFromCache(int taskId)
	{
		for (String key : assignedTasks.keySet())
			if (assignedTasks.get(key).id.equals(taskId))
				return assignedTasks.get(key);

		return taskQueue.getTaskById(taskId);
	}

	private void resubmitTask(Task task)
	{
		// an attempt to recover failure by providing more memory

		if(!isParentAlive(task.parentTaskId)){ // resubmission only works for children tasks
			removeTaskFromCache(task.id); //in case if it has not yet removed
			removeTaskFromCache(task.parentTaskId); //in case if it has not yet removed
			return;
		}

		if(task.getDetailedStatus().contains(CriticalException.CRITICAL_ERROR)) {
			task.status = "error";

			String message ="Critical error " + task + " " + task.getDetailedStatus() + "\n name:" + getOriginalName(task) + "\n" +
					task.calcServerId +  " (user " + task.getUser()+ ")";

			logger.info(message);
			Mailer.notifyDevelopers("Critical error", message);
			return;
		}

		if(task.resubmitted == null){
			task.resubmitted = 1;
			task.setMinRequiredMemory(task.getMinRequiredMemory()< 1024? 1024: task.getMinRequiredMemory() + 512);
			task.priority -= PRIORITY_DECREMENT_FOR_RESTART;
		}else{
			task.setMinRequiredMemory(task.getMinRequiredMemory() * 2);// providing dramatic increase of memory
			task.resubmitted++;
			task.priority -= PRIORITY_DECREMENT_FOR_RESTART;
		}

		/*if(task.parentTaskId != null) // update all tasks from the parent, to avoid multiple updates
			session().createQuery("update Task set minRequiredMemory = :rm, resubmitted = 1, priority = :pr where parent_task_id = :id and status = :st and minRequiredMemory < :rm").
			setInteger("id", task.parentTaskId).
			setDouble("pr", task.priority).
			setString("st", "init").
			setInteger("rm", task.getMinRequiredMemory()).
			executeUpdate();*/

		long howlong = task.timeAssigned != null ? (Calendar.getInstance().getTimeInMillis() - task.timeAssigned.getTime())/1000:0;

		String message ="Resubmitting task " + task + ", because something went wrong with server: \n" + task.getDetailedStatus() + "\n name:" + getOriginalName(task) + "\n" +
				task.calcServerId + " trial=" + task.resubmitted + " minmemory="+task.getMinRequiredMemory() + " (user "
				+ task.getUser()+ ") was calculated for " + howlong + " sec.";

		logger.info(message);
		Mailer.notifyDevelopers("Resubmitting: Server did not return task result", message);

		task.status = Task.INIT;
		task.setDetailedStatus("Resubmitted! Waiting for a free server");
		assignedTasks.remove(task.calcServerId);
		task.calcServerId = null;
		task.timeAssigned = null;
		task.timeCompleted = null;

		session()
		.createQuery(
				"update Task set status = :status, calcServerId = :csi, detailedStatus = :dstatus, timeAssigned = :ta, timeCompleted = :tc, resubmitted = :rs, minRequiredMemory = :rm where id = :id")
		.setString("status", task.status).setString("csi", task.calcServerId).setString("dstatus", null)
		.setTimestamp("ta", null).setTimestamp("tc", null).setInteger("rs", task.resubmitted).
		setInteger("rm", task.getMinRequiredMemory()).
		setInteger("id", task.id).executeUpdate();

		killChildrenTasks(task.id);
		removeTaskFromCache(task.id);

		// Update the cache
		if (taskQueue.shouldPutTask(task))
		{
			//Refetch task data
			Task refetchedTask = (Task) session().get(Task.class, task.id);
			taskQueue.putTask(refetchedTask);
			logger.info("Task " + refetchedTask + " has been refetched to the queue" + (refetchedTask.getDataSize() == 0 ? " WARNING: Data is null" : ""));
		}
	}

	public String getOriginalName(Task task) {
		if(task.parentTaskId == null) return task.taskName;
		Task refetchedTask = (Task) session().get(Task.class, task.parentTaskId);
		return getOriginalName(refetchedTask);
	}

	private Command registerServer(Command command)
	{
		if (terminatorMode)
			return new Command(Command.MS_TERMINATE, null);

		Task task = null;
		OnlinePeer peer = onlinePeers.get(command.senderId);
		peer.currentTask = null;

		if (peer.restartRequested)
		{
			peer.status = "restarting...";
			peer.restartRequested = false;
			return new Command(Command.MS_RESTART, null);
		}

		if (peer.logsRequested)
		{
			peer.logsRequested = false;
			return new Command(Command.MS_GET_LOGS, null);
		}

		// Send a command to the CS
		if (peer.commandToSend != null)
		{
			String cmd = peer.commandToSend;
			peer.commandToSend = null;
			return new Command(Command.MS_COMMAND, cmd);
		}

		ServerInfo receivedServerInfo = (ServerInfo) command.data;
		// Check for SID conflicts
		if (peer.serverInfo != null)
		{
			if (receivedServerInfo.random != peer.serverInfo.random)
			{
				peer.conflictTime = peer.lastActive;
				logger.info("Peer " + command.senderId + " has changes session ID from " + peer.serverInfo.random + " to " + receivedServerInfo.random);
			}
			if (receivedServerInfo.configurationXml == null)
				receivedServerInfo.configurationXml = peer.serverInfo.configurationXml;
		}

		peer.serverInfo = receivedServerInfo;

		if (UpdateServer.instance.currentRelease != null)
			currentVersion = UpdateServer.instance.currentRelease.version;

		if (currentVersion != null && peer.serverInfo != null && peer.serverInfo.version != null && !peer.serverInfo.version.equals(currentVersion)
				&& !peer.serverInfo.version.equals(Command.LOCAL) && !peer.serverInfo.version.equals(Command.FIXED) )
			return new Command(Command.MS_UPDATEREQUIRED, null);

		if (!isAssignNewTasksAllowed())
			return null;

		synchronized (assignedTasks)
		{
			if ((task = assignedTasks.get(command.senderId)) != null)
			{
				// Server is asking for another task, but hasn't yet completed the previous one!

				logger.info("Server " + command.senderId + " (session ID " + peer.serverInfo.random + ") did not return task result for " + task);

				assignedTasks.remove(task.calcServerId);
				removeTaskFromCache(task.id);
				if (isToBeResubmitted(task))
					resubmitTask(task);
				else
				{
					Mailer.notifyDevelopers("Metaserver: Server did not return task result", "Server " + command.senderId + " (user "
							+ task.getUser()+ ") did not return task result for " + task);

					task.status = Task.ERROR;
					task.setDetailedStatus("Server did not return task result");					
					changeTaskStatus(task);
				}
				// we will ask it to restart, just in case (maybe it has some problems, e.g. disk space is low, etc.) and also to allow this task to be taken by another server
				return new Command(Command.MS_RESTART, null); 
			}
		}

		// Update the ping map
		long pingTime = Calendar.getInstance().getTimeInMillis();
		for (String availableTask : peer.getAvailableTaskTypes()) 
			taskTypePing.put(availableTask, pingTime);

		// Actually probe a task from the queue
		synchronized (taskQueue)
		{
			task = taskQueue.getBestTaskForPeer(peer);
		}

		//peer.supportedTaskTypes = ((String) command.data).replace(",", ", ");

		// Assign task to this server
		if (task != null && !peer.disabled)
		{
			if (!task.canBeCalculatedBy(command.senderId))
			{
				if (Calendar.getInstance().getTimeInMillis() - task.time.getTime() > EXCLUSIVE_TASKS_TIMEOUT)
				{
					// This exclusive task has been waiting too long. The preferred server is not available
					taskQueue.removeTask(task);
					session().createQuery("update Task set status = :status, detailedStatus = :dstatus where id = :id").setString("status", Task.ERROR)
					.setString("dstatus", "The preferred server " + task.getPreferredServer() + " is not available for task " + task).setInteger("id", task.id)
					.executeUpdate();
					reloadCache(true, true);
				}
				peer.status = "idle";
				return null;
			}

			if (forcedTaskBindings.get(task.taskType) != null && !forcedTaskBindings.get(task.taskType).equals(command.senderId))
			{
				// This task type is reserved for a particular server
				return null;
			}

			peer.currentTask = task;
			task.status = Task.ASSIGNED;
			task.calcServerId = command.senderId;
			task.timeAssigned = new Timestamp(Calendar.getInstance().getTimeInMillis());
			task.setDetailedStatus("Assigned to " + command.senderId);

			session().createQuery("update Task set status = :status, calcServerId = :csi, detailedStatus = :dstatus, timeAssigned = :ta where id = :id")
			.setString("status", task.status).setString("csi", task.calcServerId).setString("dstatus", task.getDetailedStatus())
			.setTimestamp("ta", task.timeAssigned).setInteger("id", task.id).executeUpdate();

			logger.info("Server " + command.senderId + "(" + peer.serverInfo.random + ") has taken " + task
					+ (task.getDataSize() == 0 ? " WARNING: Data is null" : ""));

			taskQueue.removeTask(task);

			task.clearData();
			task.clearResult();

			assignedTasks.put(command.senderId, task);

			// Load a task of this task type
			queueNextTask(task);

			peer.status = "has just taken a task!";

			Task refetchedTask = (Task) session().get(Task.class, task.id); // A weird awkward way to solve the problem of not-storing task data locally in assignedTasks
			// but still be able to send this data to the server.
			// Another way would be to clone the Task somehow and put one copy to assignedTasks without data and send the other copy.

			return new Command(Command.MS_ASSIGN_TASK, refetchedTask);
		}
		else
			peer.status = "idle";

		return null;
	}

	private void queueNextTask(Task taskToReplace) {
		long time = Calendar.getInstance().getTimeInMillis();
		Task nextTask = getTask(Task.INIT, taskToReplace.taskType, taskToReplace.getMinRequiredMemory());
		if (nextTask != null)
		{
			taskQueue.putTask(nextTask);
			logger.info("Task " + nextTask + " has been added to the queue " + (nextTask.getDataSize() == 0 ? " WARNING: Data is null" : ""));
		}
		logger.info("Queue for task type " + taskToReplace.taskType + " updated in " + (Calendar.getInstance().getTimeInMillis() - time) + "ms.");
	}

	private synchronized void reloadCache(boolean reloadAssignedTasks, boolean reloadQueuedTasks)
	{
		logger.info("Reloading metaserver cache...");
		long timer = Calendar.getInstance().getTimeInMillis();

		if (reloadQueuedTasks)
		{
			// Fetch one new task for every task type
			List<Object[]> groupList = session().createCriteria(Task.class)
					.add(Restrictions.eq("status", Task.INIT))
					.setProjection(
							Projections.projectionList()
							.add(Projections.groupProperty("taskType"))
							.add(Projections.groupProperty("minRequiredMemory"))
							).list();

			logger.info("RC, check-point 1 " + (Calendar.getInstance().getTimeInMillis() - timer) + "ms.");
			List<Task> taskList = new ArrayList<Task>();
			for (Object[] groupElement : groupList)
			{
				String taskType = (String)groupElement[0];
				Integer minMemory = ((Number)groupElement[1]).intValue();
				Task waitingTask = getTask(Task.INIT, taskType, minMemory);
				if (waitingTask != null)
					taskList.add(waitingTask);
			}

			logger.info("RC, check-point 2 " + (Calendar.getInstance().getTimeInMillis() - timer) + "ms.");
			synchronized (taskQueue)
			{
				taskQueue.clear();
				for (Task task : taskList)
				{
					session().evict(task);
					taskQueue.putTask(task);
				}
			}
		}

		logger.info("RC, check-point 3 " + (Calendar.getInstance().getTimeInMillis() - timer) + "ms.");

		if (reloadAssignedTasks)
		{
			// TODO: load KILL also
			synchronized (assignedTasks)
			{
				assignedTasks.clear();
				List<Number> taskIdList = session().createCriteria(Task.class).add(
						Restrictions.eq("status", Task.ASSIGNED)).setProjection(Projections.id())
						.list();

				for (Number taskId : taskIdList)
				{
					Task task = (Task) session().get(Task.class, taskId.intValue());
					session().evict(task);
					task.clearData();
					task.clearResult();
					assignedTasks.put(task.calcServerId, task);
				}

				taskIdList = session().createCriteria(Task.class).add(
						Restrictions.eq("status", Task.STOP)).setProjection(Projections.id())
						.list();

				for (Number taskId : taskIdList)
				{
					Task task = (Task) session().get(Task.class, taskId.intValue());
					session().evict(task);
					task.clearData();
					task.clearResult();
					assignedTasks.put(task.calcServerId, task);
				}

			}
		}

		logger.info("RC, check-point 4 " + (Calendar.getInstance().getTimeInMillis() - timer) + "ms.");
		// Update list of online peers

		timer = Calendar.getInstance().getTimeInMillis() - timer;
		logger.info("Metaserver cache has been reloaded in " + timer + "ms");
	}

	private void cleanupPeerList(boolean cleanClients)
	{
		synchronized (onlinePeers)
		{
			for (Iterator<Map.Entry<String, OnlinePeer>> iter = onlinePeers.entrySet().iterator(); iter.hasNext();)
			{
				Map.Entry<String, OnlinePeer> entry = (Map.Entry<String, OnlinePeer>) iter.next();
				OnlinePeer peer = (OnlinePeer) entry.getValue();
				if (cleanClients) // Clean client peers, not too important
				{
					if (peer.getPing() > 60000 && peer.isClient)
						iter.remove();
				} else // Clean server peers that are (hopefully) overdue / disconnected
				{
					if (peer.getPing() > 2 * METASERVER_TASK_OVERDUE && !peer.isClient)
						iter.remove();
				}
			}
		}
	}

	private void onTaskOverdue(Task task)
	{
		Mailer.notifyDevelopers("Metaserver: Task overdue", "Task " + task + " calculated by " + task.calcServerId
				+ " is overdue.\nTask has been calculated for " + niceTime(task.getCalculationTime()));
		logger.info("Task " + task + " calculated by " + task.calcServerId + " is overdue. Task has been calculated for " + niceTime(task.getCalculationTime()));

		if (task.parentTaskId != null && isToBeResubmitted(task))
			resubmitTask(task);
		else
		{
			task.status = Task.ERROR;
			task.setDetailedStatus("Task is overdue");

			session().createQuery("update Task set status = :status, detailedStatus = :dstatus where id = :id").setString("status", task.status)
			.setString("dstatus", task.getDetailedStatus()).setInteger("id", task.id).executeUpdate();

			changeTaskStatus(task);
		}
	}

	private synchronized void cleanupOverdueTasks()
	{
		if (getUptime() >= 30 * 60)
		{
			logger.info("Cleaning overdue tasks...");
			List<String> disconnectedServers = new ArrayList<String>();
			for (String key : assignedTasks.keySet())
			{
				OnlinePeer peer = onlinePeers.get(key);
				if (peer == null || peer.getPing() > METASERVER_TASK_OVERDUE)
				{
					disconnectedServers.add(key);

					if (peer != null)
						peer.currentTask = null;
				}
			}

			for (String server : disconnectedServers)
			{
				//				Task task = (Task) session().get(Task.class, assignedTasks.get(server).id);
				Task task = assignedTasks.get(server);
				onTaskOverdue(task);
			}

			logger.info("Deleting the tasks scheduled for deletion and not accessed for more than 30 minutes");
			GregorianCalendar gCal = new GregorianCalendar();
			gCal.add(Calendar.HOUR, -6);
			session().createSQLQuery("delete from Task where scheduled_kill=1 and last_access < :longAgo")
			.setDate("longAgo", gCal.getTime())
			.executeUpdate();

			logger.info("Cleaning overdue finished");
		}
	}

	/**
	 * Determining list of Ids for which tasks are still in the Task 
	 * @return
	 * @throws Exception
	 */
	@Override
	public Set<String> getValidTaskReferenceIds()
	{
		Set<String> metaserverIds = new HashSet<String>();
		synchronized (this)
		{
			Transaction tx = session().beginTransaction();
			try
			{
				List<String> metaserverIdsList = session().createSQLQuery("select reference_id from Task")
						.addScalar("reference_id", StringType.INSTANCE).list();

				for(String ref:metaserverIdsList) // we will not delete old reference style data
					if(ref.length() >= Task.REFERENCE_LENGTH) // otherwise it is old style reference: will remain intact
						for(String r : ref.split("(?<=\\G.{"+Task.REFERENCE_LENGTH+"})"))
							metaserverIds.add(r);
				tx.commit();
			} catch (Exception e)
			{
				tx.rollback();
				throw e;
			}
		}
		return metaserverIds;
	}

	private Map<String, Integer> mongoDbTasksToDelete = new HashMap<String, Integer>();

	/**
	 * Determine whether we can assign new tasks based on the count of to-delete tasks in MongoDB
	 * Based on the "iron" principle with the switch-on range of [2000, 5000]  tasks to be deleted
	 * @return
	 */
	private boolean isAssignNewTasksAllowed()
	{
		int tasksToDelete = 0;
		for (String collectionName : mongoDbTasksToDelete.keySet())
			tasksToDelete += mongoDbTasksToDelete.get(collectionName);

		if (tasksToDelete > 10000)
			if (!stopTaskAssignment)
			{
				String msg = "There are " + tasksToDelete +  " MongoDB tasks to delete. Stopping assigning new tasks";
				Mailer.postMailSafely(new Email(Task.EMAIL_OCHEM,  msg, msg).useHTML());
				logger.info(msg);
				stopTaskAssignment = true;
			}

		if (tasksToDelete < 2000)
			if (stopTaskAssignment)
			{
				String msg = "There are " + tasksToDelete +  " MongoDB tasks to delete. Resuming assigning new tasks.";
				Mailer.postMailSafely(new Email(Task.EMAIL_OCHEM,  msg, msg).useHTML());
				logger.info(msg);
				stopTaskAssignment = false;
			}

		return !stopTaskAssignment;

	}

	private Task getTask(String status, String type)
	{
		session().flush();
		List<Number> tasks = session()
				.createSQLQuery("select task_id from Task force index (StatTypePri_Index) where status=:status and task_type=:taskType order by priority desc")
				.setString("taskType", type).setString("status", status).setMaxResults(1).list();

		if (tasks.size() > 0)
			return (Task) session().get(Task.class, tasks.get(0).intValue());
		else
			return null;

	}

	private Task getTask(String status, String type, Integer minMemory)
	{
		session().flush();

		List<Number> tasks = null;
		if (minMemory == null)
			tasks = session().createSQLQuery("select task_id from Task force index (StatTypePri_Index) where status=:status and task_type=:taskType and min_required_memory is null order by priority desc").setString("taskType", type).setString("status", status).setMaxResults(1).list();
		else
			tasks = session().createSQLQuery("select task_id from Task force index (StatTypePri_Index) where status=:status and task_type=:taskType and min_required_memory=:mem order by priority desc").setString("taskType", type).setString("status", status).setInteger("mem", minMemory).setMaxResults(1).list();

		if (tasks.size() > 0)
			return (Task) session().get(Task.class, tasks.get(0).intValue());
		else
			return null;

	}

	private OnlinePeer addPeer(String sid)
	{
		synchronized (onlinePeers)
		{
			OnlinePeer peer = onlinePeers.get(sid);
			if (peer == null)
			{
				peer = new OnlinePeer();
				peer.ipAddress = ServerServlet.resolveRemoteAddr(request.get());
				logger.info("New peer: " + sid);
			}else
				logger.info("Update peer: " + sid + " time: " + (Calendar.getInstance().getTimeInMillis() - peer.lastActive) + " ms");
			peer.lastActive = Calendar.getInstance().getTimeInMillis();
			onlinePeers.put(sid, peer);
			return peer;
		}
	}

	public MetaServer(String id)
	{
		super(id);
		instance = this;
		startupTime = Calendar.getInstance().getTimeInMillis();
		onlinePeers = new HashMap<String, OnlinePeer>();
		assignedTasks = new HashMap<String, Task>();
		taskQueue = new TaskQueue();

		migrateDatabase();

		Transaction tx = session().beginTransaction();

		new HangingCommandDetector().start();
		if (!mirror)
		{
			new NoSQLCleanerThread(Task.TASK_DATABASE, NoSqlTransport.DEFAULT_COLLECTION, this).start();
			new StatisticsLogger().start();
			new OverdueCleanerThread().start();
			new OnlinePeerCleanerThread().start();
			new CacheReloadThread(3600).start();
		}

		try
		{
			doInitialCleanup();
			reloadCache(true, true);
			tx.commit();
		} catch (Exception e)
		{
			tx.rollback();
			throw new RuntimeException(e);
		}
	}

	public static MetaServer getInstance()
	{
		return instance;
	}

	private void doInitialCleanup()
	{
		logger.info("Performing the initial cleanup..");
		// Kill children of killed tasks
		List<Object[]> rows = session()
				.createSQLQuery(
						"select Child.task_id, Child.parent_task_id from Task Child inner join Task Parent on (Child.parent_task_id=Parent.task_id) where (Child.status=\"init\" or Child.status=\"assigned\") and (Parent.status=\"error\" or Parent.status=\"kill\")")
				.list();
		for (Object[] row : rows)
			logger.info("[Cleanup] Killing task " + row[0] + " as a child of " + row[1]);
		if (rows.size() > 0)
			session()
			.createSQLQuery(
					"update Task Child inner join Task Parent on (Child.parent_task_id=Parent.task_id) set Child.status=\"kill\" where (Child.status=\"init\" or Child.status=\"assigned\") and (Parent.status=\"error\" or Parent.status=\"kill\")")
			.executeUpdate();

		// delete old Test tasks or those submitted by anonymous or test users
		session().createSQLQuery("delete from Task where user like \"" +Task.TEST_USER_PREFIX + "%\" or user is NULL or user =\"test\"").executeUpdate();

		// delete old crashed tasks, except for the upper level ones
		session().createSQLQuery("delete from Task where parent_task_id is not NULL and (status =\"error\" OR status=\"kill\")").executeUpdate();

		// Delete old tasks
		GregorianCalendar gCal = new GregorianCalendar();
		gCal.add(Calendar.MONTH, -4); // changed for four month to avoid error tasks in Pending tasks (which will be deleted in 3.5 months)
		Timestamp oneMonthAgo = new Timestamp(gCal.getTimeInMillis());

		List<Number> taskIds = session().createSQLQuery("select task_id from Task where time_completed < :oneMonthAgo")
				.setParameter("oneMonthAgo", oneMonthAgo).list();
		for (Number taskId : taskIds)
			logger.info("[Cleanup] Deleting an old task " + taskId);
		if (!taskIds.isEmpty())
			session().createSQLQuery("delete from Task where time_completed < :oneMonthAgo").setParameter("oneMonthAgo", oneMonthAgo).executeUpdate();

		Number oldestTaskId = (Number) session().createSQLQuery("select min(task_id) from Task where time_completed is not null").uniqueResult();
		// Added a special index to MetaServer table for this ^ query. For it the time_completed index was used, and the task_id was then obtained by full scan of results. Sometimes - up to 5 minutes.
		// P.S. There was already such an index,  I dont know why it could disappear
		if (oldestTaskId != null)
		{
			logger.info("Deleting task with id < " + oldestTaskId);
			session().createSQLQuery(
					"delete from Task where (time_completed is null) and (status <> 'init') and (status <> 'assigned') and task_id < " + oldestTaskId)
			.executeUpdate();
		}

		logger.info("The cleanup is finished.");

	}

	public long getUptime() // in seconds
	{
		return (Calendar.getInstance().getTimeInMillis() - startupTime) / 1000;
	}

	private boolean isTaskOverdue(Task task)
	{
		if (!Task.isAssigned(task.status))
			return false;
		if (task.calcServerId == null)
			return false;

		OnlinePeer calculatingPeer = onlinePeers.get(task.calcServerId);

		if (calculatingPeer == null) {
			// If the server is not the list, task is overdue if the metaserver is running long enough
			boolean res = getUptime() > METASERVER_TASK_OVERDUE / 1000; // 60 min
			if(res)logger.info("Peer is not available for task " + task.calcServerId);
			return res;
		}

		// If the server did not contact us for more than 60 min, the task is overdue
		boolean res = calculatingPeer.getPing() > METASERVER_TASK_OVERDUE;

		if(res)logger.info("Server for task: " + task.calcServerId + " did not conect during during " + METASERVER_TASK_OVERDUE/1000 + " sec");

		return res;
	}

	public static String getLocalURL(HttpServletRequest request) {
		return   request.getScheme() + "://" +
				request.getServerName() + 
				("http".equals(request.getScheme()) && request.getServerPort() == 80 || "https".equals(request.getScheme()) && request.getServerPort() == 443 ? "" : ":" + request.getServerPort() ) +
				request.getRequestURI();
	}

	public synchronized void serveGETRequest(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		if (req.getParameter("action") != null)
		{
			SimpleController.getInstance().invoke(req, res);
			return;
		}

		if (req.getParameter("server") != null)
		{
			OnlinePeer peer = onlinePeers.get(req.getParameter("server"));
			if (req.getParameter("disable") != null)
				peer.disabled = true;
			if (req.getParameter("enable") != null)
				peer.disabled = false;
			if (req.getParameter("disabletask") != null)
				peer.disableTask(req.getParameter("disabletask"));
			if (req.getParameter("enabletask") != null)
				peer.enableTask(req.getParameter("enabletask"));
			if (req.getParameter("restart") != null)
				peer.restartRequested = true;
			if (req.getParameter("getlogs") != null)
				peer.logsRequested = true;			
			if (req.getParameter("command") != null)
				peer.commandToSend = req.getParameter("command");
			res.sendRedirect(getLocalURL(req));
			return;
		}
		terminatorMode = req.getParameter("terminator") != null;

		if (UpdateServer.instance.currentRelease != null)
			currentVersion = UpdateServer.instance.currentRelease.version;
		else
			currentVersion = "";
		Transaction tx = session().beginTransaction();
		cleanupOverdueTasks();
		try
		{
			if (req.getParameter("killtask") != null)
			{
				Task task = (Task) session().get(Task.class, Integer.valueOf(req.getParameter("killtask")));
				logger.info("[Webinterface] Killing task " + task);
				task.status = Task.KILL;
				changeTaskStatus(task);
				res.setStatus(HttpServletResponse.SC_MOVED_TEMPORARILY);
				res.setHeader("Location", getLocalURL(req));
			}

			if (req.getParameter("restarttask") != null)
			{
				Task task = (Task) session().get(Task.class, Integer.valueOf(req.getParameter("restarttask")));
				logger.info("[Webinterface] Restarting task " + task);
				task.status = Task.INIT;
				changeTaskStatus(task);
				res.setStatus(HttpServletResponse.SC_MOVED_TEMPORARILY);
				res.setHeader("Location", getLocalURL(req));
			}

			if (req.getParameter("stoptask") != null)
			{
				Task task = (Task) session().get(Task.class, Integer.valueOf(req.getParameter("stoptask")));
				logger.info("[Webinterface] Requesting task to stop " + task);
				task.status = Task.STOP;
				changeTaskStatus(task);
				res.setStatus(HttpServletResponse.SC_MOVED_TEMPORARILY);
				res.setHeader("Location", getLocalURL(req));
			}

			showStatusPageInternal(req, res);
			tx.commit();
		} catch (Exception e)
		{
			e.printStackTrace();
			tx.rollback();
			throw new RuntimeException(e);
		}
	}

	public static String niceTime(long totalSec)
	{
		long tmp = totalSec;
		long sec = tmp % 60;
		long min = Math.round(Math.floor(totalSec / 60)) % 60;
		long hours = Math.round(Math.floor(totalSec / 60 / 60)) % 24;
		long days = Math.round(Math.floor(totalSec / 60 / 60 / 24));
		List<String> parts = new ArrayList<String>();
		String st = "";
		if (days > 0)
			parts.add(days + " days ");
		if (hours > 0)
			parts.add(hours + " hours ");
		if (min > 0)
			parts.add(min + " min ");
		if (sec > 0)
			parts.add(sec + " sec ");

		for (int i = 0; i < parts.size() && i < 2; i++)
			st += parts.get(i);

		return st.trim();
	}

	protected void showStatusPageInternal(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		cleanupPeerList(true);
		//		if (req.getParameter("reloadcache") != null)
		reloadCache(true, true);
		req.getSession().setAttribute(METASERVER, this);
		req.getRequestDispatcher("/jsp/metaserver-status.jsp").forward(req, res);
	}

	private void checkTaskSize(Task task)
	{
		if (task.getReferenceSize() > 900)
		{
			task.status = Task.ERROR;
			task.setDetailedStatus("Large tasks are not supported any more. Move to NoSQL storage");
		}
	}

	/**
	 * Migrate the database schema to the current version using FlyWay
	 */
	private void migrateDatabase()
	{

		Flyway flyway = Flyway.configure().dataSource(dbConfiguration.getProperty("hibernate.connection.url"), 
				dbConfiguration.getProperty("hibernate.connection.username"), dbConfiguration.getProperty("hibernate.connection.password")).load();
		flyway.baseline();
		flyway.migrate();
	}

	private class CacheReloadThread extends DaemonThread
	{
		int period;  // in seconds

		public CacheReloadThread(int seconds)
		{
			super(logger);
			this.period = seconds;
		}

		@Override
		public void wrapped() throws InterruptedException
		{
			try
			{
				Thread.sleep(period * 1000);
				Transaction e = session().beginTransaction(); // just in case we do not have an active transaction
				reloadCache(true, true); // reloading cache automatically each hour
				e.commit();
			} catch (Exception e)
			{
				logger.info(e);
				e.printStackTrace();
			}
		}
	}

	private class OverdueCleanerThread extends DaemonThread
	{
		public OverdueCleanerThread() {
			super(logger);
		}

		@Override
		public void wrapped() throws Exception
		{
			Thread.sleep(1000 * 60 * 30);

			if (MemoryUtils.getCurrentMemoryUsedFraction() > 0.9)
				Mailer.postMailSafely(new Email(Mailer.developers, "Metaserver's memory is running out", "Metaserver's memory is running out:\n" + MemoryUtils.memorySummary()).useHTML());

			Transaction tx = session().beginTransaction();
			try
			{
				cleanupOverdueTasks(); // Check for overdue tasks time to time // Midnighter on Dec 15, 2011
				tx.commit();
			} catch (Exception e)
			{
				logger.info(e);
				Mailer.notifyDevelopers(e, "Metaserver$MaintenanceThread/OverdueTasks");
				tx.rollback();
			}
		}

	}

	private class OnlinePeerCleanerThread extends DaemonThread
	{
		public OnlinePeerCleanerThread() {
			super(logger);
		}

		@Override
		public void wrapped() throws InterruptedException
		{
			Thread.sleep(1000 * 60 * 30);
			cleanupPeerList(false);
		}
	}

	// The thread periodically checks whether there is a hangind command and prints a warning
	private class HangingCommandDetector extends DaemonThread
	{
		public HangingCommandDetector() {
			super(logger);
			// TODO Auto-generated constructor stub
		}

		@Override
		public void wrapped() throws InterruptedException
		{
			Thread.sleep(2000);
			if (currentCommand != null)
			{
				long hangTime = Calendar.getInstance().getTimeInMillis() - lockStartedTime;
				if (hangTime > 1500)
					logger.info("WARNING: A hanging command " + currentCommand + " is hanging for " + (hangTime / 1000) + " seconds");
			}
		}
	}

	public Number getQueuedTasksCount()
	{
		return (Number) session().createSQLQuery("select count(*) c from Task where status='init'").addScalar("c", IntegerType.INSTANCE).uniqueResult();
	}

	public Number getPrimaryQueuedTasksCount()
	{
		return (Number) session().createSQLQuery("select count(*) c from Task where status='init' and parent_task_id is null")
				.addScalar("c", IntegerType.INSTANCE).uniqueResult();
	}

	public static String htmlAngleBrackets(String str)
	{
		str = str.replaceAll("(\\r|\\n)", "XZZZX");
		str = str.replaceAll("<config(.*?)</config>", "");
		str = str.replaceAll("<mongoDbURL>(.*?)</mongoDbURL>", "");
		str = str.replaceAll("<", "&lt;");
		str = str.replaceAll(">", "&gt;");
		str = str.replaceAll("    XZZZX", "");
		str = str.replaceAll("XZZZX", "\n");
		return str;
	}

	class StatisticsLogger extends DaemonThread
	{
		private long lastTime;
		private long lastConnectionsCount;
		private long lastLostConnectionsCount;
		private int lastNewTasks;
		private int lastCompletedTasks;
		private int lastErrors;

		@Override
		public void wrapped() throws InterruptedException
		{
			Thread.sleep(60000);
			Transaction tx = session().beginTransaction();
			try
			{
				long curTime = Calendar.getInstance().getTimeInMillis();
				int curConnectionsCount = ServerServlet.connectionAttempts.get();
				int curLostConnectionsCount = curConnectionsCount - ServerServlet.acceptedConnectionAttempts.get();
				int curNewTasks = newTasksCounter.get();
				int curCompletedTasks = completedTasksCounter.get();
				int curErrors = errorsCounter.get();

				StatisticsLog log = new StatisticsLog();
				log.date = new Timestamp(curTime);
				log.assignedTasks = Long.valueOf(assignedTasks.size());
				log.connectionsPerSecond = 1000 * (curConnectionsCount - lastConnectionsCount) / (curTime - lastTime);
				log.lostConnectionsPerSecond = 1000 * (curLostConnectionsCount - lastLostConnectionsCount) / (curTime - lastTime);
				log.newTasks = curNewTasks - lastNewTasks;
				log.completedTasks = curCompletedTasks - lastCompletedTasks;
				log.errors = curErrors - lastErrors;
				session().save(log);

				lastTime = curTime;
				lastConnectionsCount = curConnectionsCount;
				lastLostConnectionsCount = curLostConnectionsCount;
				lastNewTasks = curNewTasks;
				lastCompletedTasks = curCompletedTasks;
				lastErrors = curErrors;

				tx.commit();
			} catch (Exception e)
			{
				logger.warn("Exception in StatisticsLogger:");
				logger.warn(e);
				tx.rollback();
			}
		}

		public StatisticsLogger()
		{
			super(logger);
			lastTime = Calendar.getInstance().getTimeInMillis();
			lastConnectionsCount = ServerServlet.connectionAttempts.get();
			lastLostConnectionsCount = lastConnectionsCount - ServerServlet.acceptedConnectionAttempts.get();
			lastNewTasks = newTasksCounter.get();
			lastCompletedTasks = completedTasksCounter.get();
			lastErrors = errorsCounter.get();
		}
	}

	@Override
	public void taskToBeDeleted(String collection, int size) {
		// TODO Auto-generated method stub

	}

	@Override
	public Logger getLogger() {
		return logger;
	}

	@SuppressWarnings({ "unused" })
	public static void main(String args[]) throws Exception
	{
		MetaServer.mirror = true;
		MetaServer ms = new MetaServer(METASERVER);
		Transaction tx = ms.session().beginTransaction();

		List<Integer> ids = new ArrayList<Integer>();
		ids.add(16501);

		for(Integer i:ids) {

			List<Integer> idss = new ArrayList<Integer>();
			idss.add(i);
			Command com = new Command( null, (Serializable) i);
			Command coma = ms.clQueryTaskStatus(com, null);
			Task aa = (Task) coma.data;
			System.out.println(aa.status);
		}

		//List<ArchivedTask> tasks = (List)ms.session().createCriteria(ArchivedTask.class).list();
		//System.out.println(MemoryUtils.memorySummary());
	}

}
