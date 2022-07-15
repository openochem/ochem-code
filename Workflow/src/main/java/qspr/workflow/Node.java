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

package qspr.workflow;

import java.io.Serializable;

import javax.xml.bind.annotation.XmlTransient;

import org.apache.commons.lang.time.StopWatch;

import qspr.metaserver.ServerPool;
import qspr.metaserver.protocol.Task;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.exceptions.UserFriendlyException;

public class Node 
{
	public String id;
	public String taskName;
	public Object configuration;
	public Workflow parent;

	@XmlTransient
	public SocketSet input;

	@XmlTransient
	public SocketSet output;

	@XmlTransient
	public boolean ready = false;

	public boolean skip = false;

	@XmlTransient
	public int rank = 0;

	public void calculate() throws Exception
	{
		WorkflowNodeData wndInput = input.getWorkflowNodeData();

		if (skip)
			if (input.sockets.size() == 1 && output.sockets.size() == 1)
			{
				output.setWorkflowNodeData(wndInput); //TODO: Is this ok? We're literally reusing the same object.
				parent.out.println("WF: Skipping task "+taskName+" "+wndInput);
				ready = true;
				return;
			} else
				throw new UserFriendlyException("Task "+taskName+" with "+input.sockets.size()+" input ports and "+output.sockets.size()+" output ports selected to be skipped, but only tasks with 1 input port and 1 output port can be currently skipped.");

		StopWatch stopWatch = new StopWatch();
		stopWatch.start();
		parent.out.println("WF: Sending task "+taskName+" "+wndInput);

		Task task = new Task(taskName, (Serializable)configuration, wndInput,ServerPool.canCalculateTaskLocally(taskName)).setParentTask(parent.currentTask.get());
		configuration = null;
		input = null;
		Task responseTask = parent.myClient.calculateTask(task);

		if (responseTask.status.equals(Task.ERROR) || Task.KILLED.equals(responseTask.status))
		{
			parent.out.println("WF: Task "+responseTask+" failed");
			throw new UserFriendlyException("Task "+responseTask+" by "+responseTask.calcServerId+" failed:\n"+responseTask.getDetailedStatus());
		}

		output.setWorkflowNodeData((WorkflowNodeData)responseTask.getResult());

		ready = true;
		stopWatch.split();
		parent.out.println("WF: task "+responseTask+" calculated "+(WorkflowNodeData)responseTask.getResult() + " in " + stopWatch.toSplitString());
	}


	public boolean readyForCalculation()
	{
		boolean res = true;
		for (int sockNum = 0; sockNum < input.sockets.size(); sockNum++)
			if (input.getValue(sockNum) == null)
			{
				res = false;
				break;
			}
		return res;
	}
}

class Port
{
	Serializable content;
}
