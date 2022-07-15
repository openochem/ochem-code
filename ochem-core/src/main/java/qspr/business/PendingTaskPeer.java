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

package qspr.business;

import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.HibernateException;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.PendingTask;
import qspr.entities.PendingTask.TaskType;
import qspr.entities.User;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.protocol.Task;
import qspr.util.WrapperThread;

public class PendingTaskPeer 
{
	/**
	 * Batch update of the pending task statuses from the metaserver
	 * @return a flag that indicates whether there were any tasks to update
	 */
	public static boolean updateTaskStatuses(User user, Basket trainingSet, int limit) throws Exception
	{
		boolean anything = false;
		anything |= updateTaskStatuses(user, trainingSet, Task.ASSIGNED, limit);
		anything |= updateTaskStatuses(user, trainingSet, Task.INIT, limit);
		anything |= updateTaskStatuses(user, trainingSet, Task.STOP, limit);
		anything |= updateTaskStatuses(user, trainingSet, null, limit);

		return anything;
	}

	@SuppressWarnings("unchecked")
	public static int processMarkedTasks(boolean deleteAll) throws HibernateException, Exception {
		PendingTaskFilter filter = new PendingTaskFilter();
		filter.toBeDeleted = true;
		List<PendingTask> tasks = filter.createCriteria().list();

		for (PendingTask task : tasks)
			if (deleteAll)
				deletePendingTask(task,true);
			else
			{
				if (task.model != null)
					task.model.markAccessed();
			}

		return tasks.size();
	}

	public static void deletePendingTask(PendingTask pTask, boolean ignorePublishedTasks)
	{
		if(pTask.published){
			if(ignorePublishedTasks)return;
			throw new UserFriendlyException("This task is published. First unpublish it to delete.");
		}
		Globals.session().delete(pTask);
		if ( pTask.type == TaskType.MODEL_TRAINING && !pTask.model.recalculation)
			if (pTask.model.taskId != null)
				pTask.model.delete();
		if (pTask.type == TaskType.MODEL_TRAINING)
			pTask.model.taskId = null;
		try
		{
			getClient().killTask(pTask.taskId);
		}
		catch (Exception e)
		{
			// its ok
		}
	}

	private static CalculationClient getClient()
	{
		CalculationClient client = new CalculationClient("Pending-Tasks");
		client.setDeepSleepTime(0);

		return client;
	}

	
	/**
	 * Batch update of the pending task statuses from the metaserver
	 * @return a flag that indicates whether there were any tasks to update
	 */
	public static int countTasks(User user, Basket trainingSet, String status) throws Exception
	{
		Criteria c = Globals.session().createCriteria(PendingTask.class)
				.add(Restrictions.gt("taskId", 0));
		if (status != null)
			c.add(Restrictions.eq("status", status));
		else
			c.add(Restrictions.isNull("status"));

		c.add(Restrictions.eq("published", false));

		if (trainingSet != null)
		{
			c.createAlias("model", "m");
			c.add(Restrictions.eq("m.trainingSet", trainingSet));
		}
		c.setProjection(Projections.groupProperty("taskId"));

		return c.list().size();
	}
	
	/**
	 * Batch update of the pending task statuses from the metaserver
	 * @return a flag that indicates whether there were any tasks to update
	 */
	@SuppressWarnings("unchecked")
	public static boolean updateTaskStatuses(User user, Basket trainingSet, String status, int limit) throws Exception
	{
		long time = Calendar.getInstance().getTimeInMillis();

		Criteria c = Globals.session().createCriteria(PendingTask.class)
				.add(Restrictions.gt("taskId", 0));
		if (status != null)
			c.add(Restrictions.eq("status", status));
		else
			c.add(Restrictions.isNull("status"));
/*
		if (user != null)
		{
			c.createAlias("session", "s");
			c.add(Restrictions.eq("s.user", user));
		}
*/
		c.add(Restrictions.eq("published", false));

		if (trainingSet != null)
		{
			c.createAlias("model", "m");
			c.add(Restrictions.eq("m.trainingSet", trainingSet));
		}
		c.setProjection(Projections.groupProperty("taskId"));
		//c.setMaxResults(limit);

		List<Integer> taskIds = c.list();
		int allTasksCount = taskIds.size();
		if (taskIds.size() > limit)
		{
			Collections.shuffle(taskIds);
			taskIds = taskIds.subList(0, limit);
		}

		logger.info("Requesting " + taskIds.size() + " task statuses out of " + allTasksCount + " '"+status+"' tasks");
		CalculationClient client = new CalculationClient("Pending-Tasks");
		client.setDeepSleepTime(0);
		client.setTolerateMetaserverDown();
		List<Task> tasks = client.getTaskStatuses(taskIds);
		int i = 0;
		for (Task task : tasks)
		{
			logger.info("Updating task " + (++i) + " out of " + tasks.size());
			Globals.session().createSQLQuery("update PendingTask set status=:status, detailed_status=:detailedStatus where task_id=:taskId")
			.setParameter("status", task.status)
			.setParameter("detailedStatus", task.getDetailedStatus())
			.setParameter("taskId", task.id)
			.executeUpdate();
		}

		Globals.session().flush();
		logger.info("Task statuses for '"+status+"' updated in " + (Calendar.getInstance().getTimeInMillis() - time) + "ms.");

		return tasks.size() > 0;
	}

	public static void terminateTaskAsync(final List<Integer> taskIds) {
		if (!taskIds.isEmpty())
			new Thread(){
			@Override
			public void run() {
				logger.info("Terminating tasks " + taskIds);
				CalculationClient client = new CalculationClient("OCHEM");
				client.setTolerateMetaserverDown();

				for (Integer id : taskIds)
					try
				{
						client.killTask(id);
				} catch (Exception e)
				{
					logger.error("Could not kill the task", e);
				}
			}
		}.start();
	}

	public static void terminateTaskAsync(Integer taskId) {
		terminateTaskAsync(Arrays.asList(new Integer[]{taskId}));
	}

	public static void main(String[] args)
	{
		new WrapperThread()
		{

			@Override
			public void wrapped() throws Exception
			{
				updateTaskStatuses(null, null, 200);
			}
		}.run();
	}

	private static final transient Logger logger = LogManager.getLogger(PendingTaskPeer.class); 
}
