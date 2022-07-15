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

import java.util.ArrayList;
import java.util.List;

import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

public class SocketSet 
{
	public List<Socket> sockets;

	public SocketSet() {}

	public SocketSet(int count, boolean createSockets) 
	{
		sockets = new ArrayList<Socket>();
		for (int i = 0; i < count; i++)
			if (createSockets)
				sockets.add(new Socket());
			else
				sockets.add(null);
	}

	public DataTable getValue(int i)
	{
		return sockets.get(i).value;
	}

	public void setValue(int i, DataTable value)
	{
		sockets.get(i).value = value;
	}

	public void checkInputs() throws Exception
	{
		int i = 0;
		for (Socket socket : sockets) 
		{
			if (socket == null)
				throw new Exception("Input socket #"+i+" does not have input connector");
			i++;
		}
	}

	public WorkflowNodeData getWorkflowNodeData()
	{
		WorkflowNodeData data = new WorkflowNodeData();
		for (int i = 0; i < sockets.size(); i++)
			if (!(getValue(i) instanceof NullDataTable))
				data.setPort(i, getValue(i));
		return data;
	}

	public void setWorkflowNodeData(WorkflowNodeData data)
	{
		for (int i = 0; i < sockets.size(); i++)
			setValue(i, data.ports.get(i));	
	}
}

class Socket
{
	DataTable value;
}

// A dummy class to represent NULLs (ports, that have not been provided)
class NullDataTable extends DataTable
{
	private static final long serialVersionUID = 1L;

}
