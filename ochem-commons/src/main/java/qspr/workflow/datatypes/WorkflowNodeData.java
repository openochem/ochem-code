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

package qspr.workflow.datatypes;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import qspr.metaserver.MemoryEstimate;
import qspr.metaserver.protocol.DataSize;
import qspr.metaserver.protocol.Dimensionable;
import qspr.metaserver.protocol.Task;
import qspr.util.ClassCompressor;

import com.eadmet.exceptions.UserFriendlyException;

@XmlRootElement(name = "dataset")
public class WorkflowNodeData implements Serializable, Dimensionable, MemoryEstimate
{

	private static final long serialVersionUID = 1L;

	@XmlElement(name = "port")
	public List<DataTable> ports;

	private int nonzero = 0;

	public WorkflowNodeData()
	{
		ports = new ArrayList<DataTable>(1);
	}

	public WorkflowNodeData(DataTable dataTable)
	{
		dataTable.compact();
		ports = new ArrayList<DataTable>(1);
		ports.add(dataTable);
	}

	public WorkflowNodeData mergeTo(WorkflowNodeData arg)
	{
		if (arg != null)
		{
			for (DataTable dt : arg.ports) 
			{
				this.addPort(dt);				
			}
		}
		return this;
	}

	public WorkflowNodeData addPort(DataTable dataTable)
	{
		ports.add(dataTable);
		return this;
	}

	public void reset()
	{
		for (int i = 0; i < ports.size(); i++)
			ports.get(i).reset();
	}

	public boolean nextRow()
	{
		boolean moreRowsAny = false;
		boolean moreRowsAll = true;
		for (int i = 0; i < ports.size(); i++)
		{
			boolean moreRows = ports.get(i).nextRow();
			moreRowsAny |= moreRows;
			moreRowsAll &= moreRows;
		}

		if (moreRowsAll != moreRowsAny)
			throw new UserFriendlyException("Inconsistent WorkflowNodeData: the number of rows in different ports does not match: " + toString());

		return moreRowsAny;
	}

	public void setPort(int portNum, DataTable dataTable)
	{
		while (ports.size() <= portNum)
			ports.add(null);
		ports.set(portNum, dataTable);
	}

	public static WorkflowNodeData fromTask(Task task) throws IOException, ClassNotFoundException
	{
		return (WorkflowNodeData) task.getResult();
	}

	public DataTable getPort(String portId)
	{
		for (DataTable port : ports) {
			if (portId.equals(port.id))
				return port;
		}

		return null;
	}

	public String toString()
	{
		String res = "";
		for (int i = 0; i < ports.size(); i++)
		{
			if (ports.get(i) != null)
			{	
				String errorStr = "";

				if (ports.get(i).getRowsSize() < 10000) // A temporary workaround - do not count errors for bug datatables, save time. Later think of a more neat way
				{
					int errors = ports.get(i).errorCount();
					errorStr = (errors == 0) ? "" : "("+errors+")";
				}
				else
					errorStr = "(x)";
				String id = ports.get(i).id != null ? ports.get(i).id + ": " : "";
				res += "["+id+ports.get(i).getColumnsSize()+"x"+ports.get(i).getRowsSize()+"]"+errorStr+ " ";
			}
			else
				res += "[null] ";
		}
		return res;
	}

	public WorkflowNodeData getCopy()
	{
		// Shallow copy
		WorkflowNodeData wndCopy = new WorkflowNodeData();
		for (int i = 0; i < ports.size(); i++)
			wndCopy.addPort(ports.get(i));
		return wndCopy;
	}

	public WorkflowNodeData getDeeperCopy()
	{
		// Deeper copy
		WorkflowNodeData wndCopy = new WorkflowNodeData();
		for (int i = 0; i < ports.size(); i++)
			wndCopy.addPort(ports.get(i).getCopy());
		return wndCopy;
	}

	public WorkflowNodeData getEmptyCopy()
	{
		// Copy only the structure
		WorkflowNodeData wndCopy = new WorkflowNodeData();
		for (int i = 0; i < ports.size(); i++)
			wndCopy.addPort(ports.get(i).getEmptyCopy());
		return wndCopy;
	}

	public static WorkflowNodeData fromSerializedByteArray(byte[] data) throws Exception
	{
		return (WorkflowNodeData)ClassCompressor.byteToObject(data);
	}

	public static WorkflowNodeData fromSerializedFile(String fileName) throws Exception
	{
		byte[] data;
		File f = new File(fileName);
		data = new byte[(int)f.length()];
		FileInputStream ff = new FileInputStream(f);
		ff.read(data);
		ff.close();
		return fromSerializedByteArray(data);
	}

	// Implementation of Dimensionable interface
	public int getCols() 
	{
		int cols = 0;
		for (DataTable port : ports)
			if (port.getColumnsSize() > cols)
				cols = port.getColumnsSize();
		return cols;
	}

	public int getRows() 
	{
		if (ports == null || ports.size() == 0)
			return 0;
		int rows = 0;
		for (DataTable port : ports)
			if (port.getRowsSize() > rows)
				rows = port.getRowsSize();
		return rows;
	}

	public DataSize getValuesCount() 
	{
		DataSize size = new DataSize();
		for (DataTable port : ports){
			DataSize n = port.getValuesCount();
			size.all += n.all;
			size.nonZero += n.nonZero;
		}

		if(nonzero<size.nonZero)nonzero = (int) size.nonZero; // Not expected that this value will be long soon

		return size;
	}

	@Override
	public int getMinRequiredMemory(Task t) 
	{
		return MemoryEstimator.getMinRequiredMemory(t, this);
	}

	@Override
	public Integer getMaxNonZero() {
		if(nonzero == 0)getValuesCount();
		return nonzero;
	}	
}
