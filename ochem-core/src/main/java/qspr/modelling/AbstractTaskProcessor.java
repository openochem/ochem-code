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

package qspr.modelling;

import java.io.Serializable;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.Transactional;
import qspr.dao.Repository;
import qspr.entities.Model;
import qspr.entities.PendingTask;
import qspr.entities.PendingTask.TaskType;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.DebugLevel;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.util.WrapperThread;
import qspr.workflow.utils.QSPRConstants;

/**
 * The processor of an abstract calculation task.
 * Manages starting, registering a pending task and retrieving the ready task.
 * 
 * This class is a new concept (introduced on 24 Oct 2012), its is only implemented by ToxAlerts, SetCompare, MolOptimiser (June 2013).
 * In future, we would need to move ModelProcessor under the same hood.
 * 
 * @author midnighter
 *
 */
abstract public class AbstractTaskProcessor extends WrapperThread
{

	private static transient final Logger logger = LogManager.getLogger(AbstractTaskProcessor.class);

	public TaskType taskClass;
	public int defaultTaskPriority = TaskPriority.NORMAL;

	/**
	 * The preferred calculation server for the submitted task
	 */
	public String preferredServer = OCHEMConfiguration.defaultPreferredServer;

	public boolean exitAfterPosting = false;

	/**
	 * The description of the dataset
	 */
	public String setDescription = "";

	public String taskName = null;

	public int taskDebugLevel;

	/**
	 * A flag for saving the pending task early on, before the actual calculation task is prepared and sent out to metaserver.
	 */
	public boolean savePendingTaskEarly = false;

	private CalculationClient client;

	protected boolean taskReceived = false;

	/**
	 * Tolerate the task failures. 
	 * This is useful for cache-aware task types, 
	 * which should not fail if there is cache and all non-cached molecules failed
	 */
	public boolean allowFailures;

	/**
	 * Provides a description of the task, i.e., used conditions
	 */
	public String taskDescription;

	public AbstractTaskProcessor()
	{
		updateSessionTime = false;
	}

	/**
	 * The simplistic execution scenario
	 */
	@Override
	final public void wrapped() throws Exception
	{
		if (pTask != null)
			restoreFromPendingTask(pTask);

		runTask();

		if (!exitAfterPosting && !isError())
			setStatus("Finished");
	}

	/**
	 * Create a calculation task ready for posting to the metaserver
	 * Apparently, this is a remote task
	 * @return
	 * @throws Exception
	 */
	public Task createTask() throws Exception
	{
		return new Task(getTaskType(), getTaskConfiguration(), getTaskData(), false);
	}

	/*
	 * A model related to the pending task. Null by default.
	 */
	protected Model getModel()
	{
		return null;
	}

	/**
	 * Required for synchronisation block
	 */

	private static Integer taskSynchronise = 0, allTask = 0;

	/**
	 * Start a calculation task and track its status
	 * @throws Exception
	 */
	public void runTask() throws Exception
	{
		// Register a pending task early with taskId = 0


		if (savePendingTaskEarly)  // creates a stub of pendingTask
		{
			int priority = TaskPriorityManager.getNewTaskPriority(getTaskPriority(), pTask, Globals.userSession());		
			if (pTask == null){
				pTask = new PendingTask(taskClass, 0, getAttachment()).setPriority(priority).setDescription(getSetDescription()).setModel(getModel());
				logger.info("Pending task " + taskClass +" is created without an attachment...");
			}

			pTask.taskId = 0; // PseudoId
			pTask.name = taskName;
			pTask.setPriority(priority);
			registerPendingTask(pTask);
			logger.info("Pending task " + taskClass +"  is registered early...");
		}

		setStatus("Waiting to start preparation of a calculation task...");

		Globals.commitAllTransactions(); // for a moment we will wait

		int taskID = -1;  // not yet defined

		synchronized (taskSynchronise) {
			Globals.startAllTransactions();

			setStatus("Preparing the calculation task...");
			logger.info("entering the synchronized block  ... " + allTask);

			CalculationClient client = getClient();
			client.setTolerateMetaserverDown();

			logger.info("Preparing the calculation task...");
			Task task = createTask();

			if(task !=null && task.taskName == null)
				task.taskName = taskName;

			int priority = 0;

			if(!Globals.isStandalone){ 
				if(task != null && task.datarows != null)
					Repository.user.checkEligibility(task.datarows, bonusMultiplier());

				priority = TaskPriorityManager.getNewTaskPriority(getTaskPriority(), pTask, Globals.userSession());
			}

			// Task can be null, then its a dummy pending task with no real calculations
			if (task != null)
			{
				task.setPriority(priority);
				task.debug = taskDebugLevel;
				if (preferredServer != null && !"".equals(preferredServer.trim())) {
					task.setPreferredServer(preferredServer);
					if(task.debug == DebugLevel.NONE)
						task.debug = DebugLevel.ALL;
				}

				logger.info("Posting the task");
				taskID = client.postTask(task);
				logger.info("Task has been posted with id "+taskID);
			}

			if (pTask == null){  // Used for ToxAlerts screening
				//TODO this part of code in particular getAttachment() significantly slow down use f FeatureNet tasks: all molecules are stored in the database!
				pTask = new PendingTask(taskClass, taskID, getAttachment()).setPriority(priority).setDescription(getSetDescription()).setModel(getModel());
				pTask.name = taskName;
				logger.warn("Second path to create models is used.");
			}
			else
			{
				if(pTask.attachment == null)
					pTask.setAttachment(getAttachment()); // model does not have yet attachment...
				pTask.taskId = taskID;
				pTask.description = taskDescription;
			}

			registerPendingTask(pTask); // saving the state to the database

			logger.info("Finishing the synchronized block  ... " + allTask++);

			Globals.commitAllTransactions();
		}

		onTaskPosted();

		logger.info("Task has being posted ...");

		if (taskID != -1 && !exitAfterPosting)
		{
			while (!update())
				Thread.sleep(1000);
			logger.info("Task is finished  ...");
		}
	}

	public int bonusMultiplier()
	{
		return QSPRConstants.MODEL_BONUS;
	}

	/**
	 * This method is mostly intended for cleaning up memory after the task has been posted.
	 * Override it if necessary
	 */
	public void onTaskPosted() {
		// Nothing to clean up by default
	}

	public boolean isError()
	{
		return getStatus() != null && getStatus().toLowerCase().startsWith("error");
	}

	/**
	 * Update the status of the processor. If required, fetch the status of a task from metaserver
	 */
	public synchronized boolean update() throws Exception
	{
		if (isReady())
			return true;

		if (pTask != null && pTask.taskId != 0)
		{
			CalculationClient client = getClient();

			final Task task = client.getTask(pTask.taskId);

			if (task != null)
			{
				try
				{
					if (!allowFailures)
						task.check();

					new Transactional() {
						@Override
						protected void wrapped() throws Exception {
							setStatus("Task received. Preparing to display results...");
							onTaskReceived(task);
						}
					}.execute();

					setStatus(client.getStatus());
				}
				catch (Exception e)
				{
					setStatus("Error: " + e.getMessage());
					e.printStackTrace();
				}
				finally
				{
					taskReceived = true;
				}
				return true;
			} else
			{
				Task status = client.getTaskStatus(pTask.taskId, false);
				pTask.update(status);
			}
			setStatus(client.getStatus());
		}
		return false;
	}

	/**
	 * The task is ready if one of the following:
	 *  - the calculation task has been received form metaserver
	 *  - the thread is still running
	 *  - it was a dummy calculation task (e.g., all results cached)
	 */
	public boolean isReady()
	{
		// We haven't posted the real task yet. Check if the thread is running
		if (pTask == null || pTask.taskId == 0)
		{
			return !isRunning();
		}

		// Its not a real calculation task. Its already done.
		if (pTask.taskId == -1)
			return true;

		return taskReceived;
	}

	public void restoreFromPendingTask(PendingTask pTask) throws Exception
	{
		setStatus("Preparing the data");
		this.pTask = pTask;

		setDescription = pTask.setDescription;
		restoreFromAttachment(pTask.attachment.getObject());

		// Try to fetch the task result
		if (pTask.taskId != null && pTask.taskId != 0)
		{
			Task task = pTask.retrieveTask(false);
			if (task.isReady() && !task.isError())
				onTaskReceived(task);
		}
	}

	/**
	 * Defines the task priority. Can be overridden to provide customized priority restrictions
	 * @return
	 */
	public int getTaskPriority()
	{
		return defaultTaskPriority;
	}

	protected String getSetDescription()
	{
		return setDescription;
	}

	private CalculationClient getClient()
	{
		if (client == null)
		{
			client = new CalculationClient(taskClass + "/" + Globals.getClientSID(), Globals.getUsername());
			client.setTolerateMetaserverDown();
			client.setDeepSleepTime(1000);
		}
		return client;
	}

	/**
	 * The metaserver-related task type
	 * @return
	 */
	abstract protected String getTaskType();

	/**
	 * The configuration for the posted task
	 * @return
	 * @throws Exception
	 */
	abstract protected Serializable getTaskConfiguration() throws Exception;

	abstract protected Serializable getTaskData() throws Exception;

	/**
	 * The object attached to the pending task that should contain all required info to restore the state
	 * @return
	 * @throws Exception 
	 */
	abstract protected Serializable getAttachment() throws Exception;

	/**
	 * When the task is received, this method is called
	 * @param task
	 * @throws Exception
	 */
	abstract protected void onTaskReceived(Task task) throws Exception;

	/**
	 * Restore the state of the processor from the object attached to the pending task
	 * @param attachment
	 */
	abstract protected void restoreFromAttachment(Serializable attachment);

}
