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

package qspr.workflow.utils;

import java.io.IOException;
import java.io.Serializable;

import qspr.exceptions.CalculationException;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.Event;
import qspr.metaserver.EventListener;
import qspr.metaserver.ServerPool;
import qspr.metaserver.cs.WorkflowNodeServer;
import qspr.metaserver.protocol.Task;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.exceptions.AssertionFailure;
import com.eadmet.utils.MemoryUtils;

public class CalculationTask
{

	public final static int RETRIEVE_ATTEMPTS = 1000;

	// data that were sent for calculation
	public WorkflowNodeData wndInput;
	// results of calculations
	protected WorkflowNodeData wndOutput;
	public int taskId;
	protected boolean ready = false;
	public Task resultTask;
	public Object configuration;
	public String taskName;
	public boolean lazyTaskRetrieval = false;
	public String status;
	private String id = "Bag";

	protected WorkflowNodeServer parent;
	protected CalculationClient client = new CalculationClient();
	public boolean allowFailures = false;
	public String error;
	public int size;

	public boolean isReady() throws Exception
	{
		if (ready)
			return true;
		if (wndOutput != null || resultTask != null)
			return true;
		Task task = client.getTaskStatus(taskId, !lazyTaskRetrieval);
		if (task != null)
			status = task.getDetailedStatus();
		if (task != null && task.isReady())
		{
			setResultTask(task);
			parent.out.println(id + " ("+task+") finished by " + resultTask.calcServerId + ": " 
					+ MemoryUtils.memorySummary());
			onReady();
			return ready = true;
		}

		return false;
	}
	
	public void cleanup()
	{
		
	}

	private void setResultTask(Task resultTask) throws CalculationException, ClassNotFoundException, IOException
	{
		this.resultTask = resultTask;
		if (!allowFailures)
			resultTask.check();
		if (resultTask.isError())
			error = resultTask.getDetailedStatus();

		if (resultTask.hasResult())
			wndOutput = WorkflowNodeData.fromTask(resultTask);

		resultTask.clearAll();
	}

	public WorkflowNodeData getWndOutput() throws IOException, ClassNotFoundException, CalculationException
	{
		if (ready & wndOutput == null)
		{
			Task task = null;
			// Safe-fetching with multiple fetch attempts should be implemented on the level of CalculationClient.
			// See CalculationClient.postTask for example. / Midnighter on Jun 6, 2011
			// This machinery should be transparent and invisible form here.
			for (int a = 0; a < RETRIEVE_ATTEMPTS; a++)
			{
				parent.out.println("Lazy fetching a task result for " + taskId + (a == 0 ? "" : (" attempt=" + a)));
				task = client.getTask(taskId);
				if (task != null)
					break; // OK! retrieved!
			}

			if (task == null)
				throw new IOException("Failed to retrieve task results after " + RETRIEVE_ATTEMPTS + " attempts");

			setResultTask(task);

			if (wndOutput == null)
				throw new AssertionFailure("Task output is null");
		}

		return wndOutput;
	}

	public void clean()
	{
		wndOutput = null;
	}

	public void post() throws IOException, ClassNotFoundException, InterruptedException{
		post(true);
	}

	
	@SuppressWarnings({ "unchecked", "rawtypes" })
	public void post(boolean firstCalculateLocally) throws IOException, ClassNotFoundException, InterruptedException
	{
		size = wndInput.ports.get(0).getRowsSize();
		client.setTolerateMetaserverDown();
		client.statusChange.addListener(new EventListener()
		{
			public void onEvent(Event event, Object arg)
			{
				if (client.getStatus() != null)
					parent.setStatus(id + ": " + client.getStatus());
			}
		});

		parent.out.println("Posting task " + id + ": " + wndInput);
		taskId = client.postTask(new Task(taskName, (Serializable) configuration, wndInput, firstCalculateLocally && ServerPool.canCalculateTaskLocally(taskName)) //DEBUG should be faster, no serialisation!
				.setParentTask(parent.currentTask.get()), firstCalculateLocally);
		configuration = null;

	}

	public void kill() throws IOException, ClassNotFoundException
	{
		client.killTask(taskId);
	}

	public CalculationTask setParent(WorkflowNodeServer parent)
	{
		this.parent = parent;
		return this;
	}

	public void onReady() throws Exception
	{

	}

	public CalculationTask setId(String id)
	{
		this.id = id;
		return this;
	}

}
