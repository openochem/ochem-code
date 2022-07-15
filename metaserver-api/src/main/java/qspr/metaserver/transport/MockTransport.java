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

import java.io.IOException;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.Task;

/**
 * A simple mock transport that returns a predefined task result.
 * Used by the tests that need to easily mock metaserver connections
 * 
 * @author midnighter
 */
public class MockTransport implements Transport, Serializable
{
	/**
	 * A predefined task result
	 */
	public Serializable defaultResult;
	
	public Task inTask;
	
	public boolean alwaysReady = false;
	
	/**
	 * Force task failure
	 */
	public boolean fail = false;
	
	public MockTransport(Serializable defaultResult)
	{
		this.defaultResult = defaultResult;
	}
	
	public MockTransport()
	{
		
	}
	
	public Serializable getTask(Task inTask) throws IOException, ClassNotFoundException
	{
		return defaultResult;
	}
	
	public Serializable getTaskResult(Task inTask) throws IOException, ClassNotFoundException
	{
		return defaultResult;
	}
	
	@SuppressWarnings("rawtypes")
	@Override
	public Command executeCommand(Command command) throws MalformedURLException, IOException, ClassNotFoundException
	{
		logger.info("Calling a transport command " + command);
		Task task;
		switch (command.id)
		{
			case Command.CL_SUBMIT_TASK:
				task = (Task) command.data;
				this.inTask = task;
				task.id = 10000000 + Integer.valueOf((int)Math.round(Math.random()*1000000));
				task.calcServerId = "MockServer";
				task.status = Task.INIT;
				return new Command(Command.MS_TASK_STATUS, task);
			case Command.CL_QUERY_TASK:
			case Command.CL_QUERY_TASK_STATUS:
				task = new Task();
				if (inTask == null && !alwaysReady)
					task.status = "init";
				else
				{
					task.timeCompleted = task.timeAssigned = new Timestamp(Calendar.getInstance().getTimeInMillis());
					task.setResult(getTaskResult(inTask));
					task.status = fail ? "error" : "ready";
					if (fail)
						task.setDetailedStatus("Simulated failure of a calculation task");
				}
				
				if (command.data instanceof List)
				{
					ArrayList<Task> tasks = new ArrayList<Task>();
					for (int i = 0; i < ((List) command.data).size(); i++)
						tasks.add(task);
					return new Command(Command.MS_TASK_STATUS, tasks);
				}
				return new Command(Command.MS_TASK_STATUS, task);
			default:
				return null;
		}
	}
	
	private static final long serialVersionUID = 1L;
	private static final Logger logger = LogManager.getLogger(MockTransport.class);
}
