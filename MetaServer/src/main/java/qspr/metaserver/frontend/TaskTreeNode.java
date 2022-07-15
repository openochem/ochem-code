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

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;

import qspr.metaserver.ArchivedTask;
import qspr.metaserver.MetaServer;
import qspr.metaserver.protocol.Task;

/**
 * A node in a tree of tasks
 * @author midnighter
 *
 */
public class TaskTreeNode
{
	/**
	 * The task of this node
	 */
	public Task task;

	/**
	 * Total time of the task and all its children
	 */
	public long totalTime;

	/**
	 * Maximal columns that has ever participated in the task or its children
	 */
	public long maxColumns;

	public List<TaskTreeNode> children = new ArrayList<TaskTreeNode>();

	public TaskTreeNode()
	{

	}

	public TaskTreeNode(Task task)
	{
		this.task = task;
		getChildren();
	}

	public TaskTreeNode(ArchivedTask task)
	{
		this.task = getArchivedTask(task);
		getChildren();
	}

	@SuppressWarnings("unchecked")
	public void getChildren()
	{
		totalTime = task.getCalculationTime();
		maxColumns = task.cols  == null ? 0: task.cols;

		List<Task> tasks = MetaServer.getInstance().session().createCriteria(Task.class).add(Restrictions.eq("parentTaskId", task.id)).addOrder(Order.desc("id")).setMaxResults(MetaServer.LIMITCHILDREN).list();

		Set<Integer> taskIDs = new HashSet<Integer>();
		for (Task childTask : tasks)
		{
			TaskTreeNode node = new TaskTreeNode();
			node.task = childTask;
			taskIDs.add(childTask.id);
			children.add(node);
		}

		if(children.size()<MetaServer.LIMITCHILDREN) {
			List<ArchivedTask> archivedTasks = MetaServer.getInstance().session().createCriteria(ArchivedTask.class).add(Restrictions.eq("parentTaskId", task.id)).addOrder(Order.desc("id")).setMaxResults(MetaServer.LIMITCHILDREN - children.size()).list();
			for (ArchivedTask childTask : archivedTasks)
			{
				if (!taskIDs.contains(childTask.id))
				{
					TaskTreeNode node = new TaskTreeNode();
					node.task = getArchivedTask(childTask);
					children.add(node);
				}
			}
		}

		for (TaskTreeNode child : children)
		{
			child.getChildren();
			totalTime += child.totalTime;
			maxColumns = Math.max(maxColumns, child.task.cols == null?0:child.task.cols);
		}
	}

	private Task getArchivedTask(ArchivedTask aTask)
	{
		Task task = new Task();
		task.id = aTask.id;
		task.parentTaskId = aTask.parentTaskId;
		task.status = aTask.status;
		task.client = aTask.client;
		task.datarows = aTask.datarows;
		task.cols = aTask.cols;
		task.setDetailedStatus(aTask.getDetailedStatus());
		task.taskType = aTask.taskType;
		task.taskName = aTask.taskName;
		task.calcServerId = aTask.calcServerId;
		task.time = aTask.time;
		task.timeAssigned = aTask.timeAssigned;
		task.timeCompleted = aTask.timeCompleted;
		task.priority = aTask.priority;
		task.referenceId = aTask.referenceId;
		task.setUser(aTask.getUser());
		task.peakMemoryUsage = aTask.peakMemoryUsage;
		task.md5 = aTask.md5;
		return task;
	}
}
