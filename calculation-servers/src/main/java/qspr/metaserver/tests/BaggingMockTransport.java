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

package qspr.metaserver.tests;
import java.io.IOException;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.Transport;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

public class BaggingMockTransport implements Transport, Serializable
{

	WorkflowNodeData defaultResult;

	public BaggingMockTransport(WorkflowNodeData defaultResult)
	{
		this.defaultResult = defaultResult;
	}

	Map<Integer, Integer> sizeMap = new HashMap<Integer, Integer>();

	public Serializable getTaskResult(Integer taskId) throws IOException, ClassNotFoundException
	{
		int trimSize = sizeMap.get(taskId);
		WorkflowNodeData wnd = new WorkflowNodeData();
		for (int i=0; i<defaultResult.ports.size(); i++)
		{
			DataTable dt = defaultResult.ports.get(i);
			DataTable newDt = dt.getCopy();
			while (newDt.getRowsSize() > trimSize)
				newDt.deleteRow(newDt.getRowsSize() - 1);
			wnd.addPort(newDt);
		}
		return wnd;
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
			WorkflowNodeData taskData = (WorkflowNodeData)task.getData();
			task.id = 10000000 + Integer.valueOf((int)Math.round(Math.random()*1000000));
			task.calcServerId = "MockServer";
			task.status = Task.INIT;
			sizeMap.put(task.id, taskData.ports.get(0).getRowsSize());
			return new Command(Command.MS_TASK_STATUS, task);
		case Command.CL_QUERY_TASK:
		case Command.CL_QUERY_TASK_STATUS:
			task = new Task();
			Integer taskId = (Integer)command.data;
			task.timeCompleted = task.timeAssigned = new Timestamp(Calendar.getInstance().getTimeInMillis());
			task.setResult(getTaskResult(taskId));
			task.status = "ready";

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
	private static final Logger logger = LogManager.getLogger(BaggingMockTransport.class);
}
