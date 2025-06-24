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

import java.util.ArrayList;
import java.util.List;
import java.util.Vector;

import qspr.exceptions.CalculationException;
import qspr.metaserver.cs.WorkflowNodeServer;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.transport.CSTransport;
import qspr.metaserver.transport.TransportFactory;

import com.eadmet.utils.MemoryUtils;

/**
 *  Provides calculation of a set of tasks (any number of more than 2 tasks)
 *  Avoid creation of an overload of metaserver 
 * @author Midnighter/itetko
 *
 */

public class CalculationTaskSet 
{
	public static float TOLERATE_ERRORS = 10f; // we will try to fetch task 10% more the number of tasks -- should work!
	public String prefix = "";
	private WorkflowNodeServer server;
	public List<CalculationTask> tasks = new ArrayList<CalculationTask>();
	public boolean allowFailedTasks = false;

	public CalculationTaskSet(WorkflowNodeServer server)
	{
		this.server = server;
	}

	/**
	 * @param zeroTrainData  -- delete train data once task is finished to decrease memory usage 
	 * @throws Exception 
	 */
	public void calculate(boolean zeroTrainData) throws Exception{
		calculate(zeroTrainData,0);
	}

	/**
	 * @param zeroTrainData  -- delete train data once task is finished to decrease memory usage 
	 * @param overrun -- stop even overrun tasks are still running
	 * @throws Exception 
	 */
	public void calculate(boolean zeroTrainData, int overrun) throws Exception 
	{
		Vector<Integer> notfinishedTasks = new Vector<Integer>();

		for (int i = 0; i < tasks.size(); i++)
		{
			notfinishedTasks.add(i);
			if (zeroTrainData)
				tasks.get(i).wndInput = null;
		}

		String message = "";

		while (notfinishedTasks.size() > overrun)
		{
			try
			{
				if (tasks.size() > 1)
					server.setStatus(prefix + "At least " + (tasks.size() - notfinishedTasks.size()) + " of " + tasks.size() + " tasks are done. " + message);
				else if (tasks.size() == 1)
					if (tasks.get(0).status != null)
						server.setStatus(prefix + tasks.get(0).status);

				int waitingTime = CSTransport.metaserverDown ? 6 * TransportFactory.waitingTimeBetweenRequests : TransportFactory.waitingTimeBetweenRequests;

				if (tasks.size() == 1 && tasks.get(0).taskId < 0)
					waitingTime = 1000;

				Thread.sleep(waitingTime);

				// This loop is required to check task faster and will lower load of the server, particular in apply model.

				for ( int n = 0; notfinishedTasks.size() > overrun && n < notfinishedTasks.size(); n++)
				{ 
					boolean manyLeft = notfinishedTasks.size() > overrun * 2;
					int i = manyLeft  ? (int) (Math.random()*notfinishedTasks.size()) : n;
					int taskNum = notfinishedTasks.get(i);

					if (!tasks.get(taskNum).isReady()){
						message = tasks.get(taskNum).status;
						message = (message != null && message.length() > 0 && !message.contains(Command.METASERVER_DOWN))? "Exemplary message: "+message : "";
						if(manyLeft) break;
						continue; // we are near to the end; let us check all the remaining tasks
					}

					if (tasks.get(taskNum).resultTask != null)
						if (!allowFailedTasks)
							tasks.get(taskNum).resultTask.check();
					tasks.get(taskNum).cleanup();
					server.out.println(MemoryUtils.memorySummary());
					notfinishedTasks.remove(i);
					n--; // start from the same number, since we have decreased the size
				}
			}
			catch (Exception e)
			{
				server.out.print("Calculation task set failed: ");
				e.printStackTrace(server.out);
				// We failed. Kill all tasks
				for (CalculationTask task : tasks) 
				{
					if (task.resultTask == null)
						task.kill();
				}
				throw e;
			}
		}
		checkFailures();
		server.setStatus("All " + tasks.size()+ " tasks have been finished");
	}

	/**
	 * Check if there is at least one successfully finished task
	 */
	private void checkFailures() throws CalculationException {
		boolean anythingSucceded = false;
		for (CalculationTask cTask : tasks)
			anythingSucceded |= (cTask.error == null);

		if (!anythingSucceded)
			throw new CalculationException(tasks.get(0).error);
	}

	public void addTask(CalculationTask cTask)
	{
		cTask.parent = server;
		tasks.add(cTask);
	}

	public CalculationTask getTask(int n)
	{
		return tasks.get(n);
	}


	public void post() throws Exception
	{
		server.out.println(prefix + "Posting " + tasks.size() + " tasks");
		// Post all tasks
		for (CalculationTask calculationTask : tasks)
		{
			calculationTask.allowFailures = allowFailedTasks;
			calculationTask.post();
		}
	}
}
