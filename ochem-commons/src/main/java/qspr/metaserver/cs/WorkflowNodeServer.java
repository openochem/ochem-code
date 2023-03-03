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

package qspr.metaserver.cs;

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.MemoryUtils;

import qspr.configuration.JAXBContextFactory;
import qspr.dao.Various;
import qspr.metaserver.CalculationServer;
import qspr.metaserver.protocol.DataSize;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.DebugLevel;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.CalculationTask;
import qspr.workflow.utils.CalculationTaskSet;
import qspr.workflow.utils.ExcelWriter;
import qspr.workflow.utils.QSPRConstants;

/**
 * An abstract calculation server that accepts and returns a set of data tables (i.e. "ports" in WorkflowNodeData.java)
 * 
 * This server supports the concept of "flows".
 * A flow is a relationship which means row-to-row correspondence between two ports.
 * That means that:
 * a) the number of rows should be same
 * b) error rows can be safely removed from input and added as errors after calculation
 * 
 * @author midnighter
 *
 */
public abstract class WorkflowNodeServer extends CalculationServer
{
	public static final String STUBROW = "stub-row";

	private HashMap<Integer, Integer> inputFlows = new HashMap<Integer, Integer>();
	private HashMap<Integer, Integer> outputFlows = new HashMap<Integer, Integer>();
	public ThreadLocal<Task> currentTask = new ThreadLocal<Task>();
	final public static JAXBContext jContext;
	private long taskStartedTime;
	int debug = DebugLevel.NONE;

	public HashMap<String, String> params = new HashMap<String, String>();

	// Default maximum size of a locally calculated task. 
	// If the task is larger than this value, it will splitted into several smaller ones, which will be reposted in parralel
	// Zero stands for disabled parallelisation 
	public int repostSize = 0;

	static private Boolean isRunningTest = null;

	static public boolean isRunningTest() {
		if (isRunningTest == null) {
			isRunningTest = true;
			try {
				Class.forName("org.junit.Test");
			} catch (ClassNotFoundException e) {
				isRunningTest = false;
			}
		}
		System.out.println("Runing junit test: " + isRunningTest);
		return isRunningTest;
	}

	static 
	{
		try {
			jContext = JAXBContextFactory.get("qspr.workflow.datatypes:qspr.metaserver.configurations");
		} catch (JAXBException e) {
			e.printStackTrace();
			throw new ExceptionInInitializerError(e);
		}
	}

	protected int getRepostSize(Serializable configuration)
	{
		// Allow servers to adjust the reposting size depending on the configuration
		return repostSize;
	}

	@Override
	public void calculate(Task task) throws Exception
	{
		setStatus("Task started: " + task.taskType);

		WorkflowNodeData wndOut=implementCalculateTask(task);
		Runtime.getRuntime().gc();

		DataSize counts = wndOut.getValuesCount();

		for (int i = 0; i < wndOut.ports.size(); i++)
		{
			if(wndOut.ports.get(i) == null) continue;
			out.println("port "+i+" columns="+wndOut.ports.get(i).getColumnsSize()+" rows="+wndOut.ports.get(i).getRowsSize());
		}

		setStatus("Serializing table with "+counts + " objects; available memory="+MemoryUtils.getCurrentMemoryFree()+" MB");

		try{
			task.setResult(wndOut);
		} catch (Throwable t) {
			t.printStackTrace();
			if(t.getMessage() != null && t.getMessage().contains("Requested array size exceeds VM limit"))
				throw new UserFriendlyException("Calculation failed to return table with "+counts +" entries: Java VM cannot handle such large data. Decrease your datasize!");
			else
				throw new UserFriendlyException(t.getMessage());
		}

	}

	/**
	 * This step was added to implicitly release memory before compression of the results
	 * @param task
	 * @return
	 * @throws Exception
	 */

	@SuppressWarnings("unchecked")
	private WorkflowNodeData implementCalculateTask(Task task) throws Exception
	{
		Various.molecule = null; // should be determined dynamically, if required
		WorkflowNodeData wndOut = null;

		try {

			taskStartedTime = Calendar.getInstance().getTimeInMillis();
			resetMaxUsedMemory();
			currentTask.set(task);
			WorkflowNodeData in = (WorkflowNodeData)task.getData();
			out.println(": -> "+in);
			checkPortsSize(new HashMap<Integer,Integer>(), in.ports, inputFlows);

			HashMap<Integer, Integer> inputFlows = (HashMap<Integer, Integer>) this.inputFlows.clone();
			for (int i = 0; i < in.ports.size(); i++)
				if (in.ports.get(i) == null)
					inputFlows.remove(Integer.valueOf(i));

			// Define common error rows for each flow group
			Map<Integer, List<ErrorRow>> errorIndices = new HashMap<Integer, List<ErrorRow>>();

			DataTable dtin = in.ports.get(0);

			if(dtin.columns.contains(QSPRConstants.SDF_COLUMN))
				Various.molecule = Various.getCheminfImpl(dtin.getChemInfEngine());

			Integer group;
			for (int input = 0; input < in.ports.size(); input++)
				if ((group = inputFlows.get(input)) != null && in.ports.get(input) != null)
				{
					List<ErrorRow> indices = errorIndices.get(group);
					if (indices == null)
					{
						indices = new ArrayList<ErrorRow>();
						errorIndices.put(group, indices);
					}

					DataTable dt = in.ports.get(input);

					int index;
					dt.reset();
					while (dt.nextRow())
						if (Task.ERROR.equals(dt.getCurrentRow().status))
						{
							if ((index = indices.indexOf(new ErrorRow(dt.currentRow))) == -1)
								indices.add(new ErrorRow(dt.currentRow).setOriginalRow(dt.getCurrentRow()));
							else
								indices.get(index).setOriginalRow(dt.getCurrentRow());
						}
					//else if (WorkflowNodeServer.STUBROW.equals(dt.getCurrentRow().status)

				}

			// Delete error rows from datatables
			List<DataTable> processedInputs = new ArrayList<DataTable>();
			for (int input = 0; input < in.ports.size(); input++)
				if ((group = inputFlows.get(input)) != null)
				{

					// Check if we already processed this datatable
					if (processedInputs.contains(in.ports.get(input)))
						continue;

					List<ErrorRow> indices = errorIndices.get(group);
					Collections.sort(indices);
					for (int i = indices.size() - 1; i >= 0; i--)
					{
						ErrorRow er = indices.get(i);
						er.setOriginalRow(in.ports.get(input).getRow(er.index));
						in.ports.get(input).deleteRow(er.index);
					}

					processedInputs.add(in.ports.get(input));
				}

			// - Stubs
			List<List<Integer>> stubs = new ArrayList<List<Integer>>();
			for (int i = 0; i < in.ports.size(); i++)
				if (inputFlows.get(i) != null)
					stubs.add(removeStubs(in.ports.get(i)));

			this.out.println("After flow correction -> "+in);
			Runtime.getRuntime().gc();

			if (currentTask.get().debug > DebugLevel.NONE)
			{
				debug = currentTask.get().debug;
				Marshaller marshaller = jContext.createMarshaller();
				marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
				marshaller.marshal(in, new File(getDebugDirectory(currentTask.get().id) + "/in-" + currentTask.get().id + "-" + supportedTaskType + ".xml"));
				new ExcelWriter(getDebugDirectory(currentTask.get().id) + "/in-" + currentTask.get().id + "-" + supportedTaskType + ".xls").addWnd(in).finalize();
			}

			Serializable configuration = task.getConfiguration();
			int suggestedRepostSize = getRepostSize(configuration);
			if (suggestedRepostSize == 0 || in.ports.size() > 1 || (task.numberCrashes() == 0 && in.ports.get(0).getRowsSize() <= suggestedRepostSize) )
				wndOut = calculate(in.getDeeperCopy(), configuration);
			else
				wndOut = calculateInParralel(in, configuration, task.numberCrashes());

			Runtime.getRuntime().gc();

			// + Stubs
			for (int i = 0; i < in.ports.size(); i++)
				if (inputFlows.get(i) != null)
					restoreStubs(in.ports.get(i), stubs.get(i));

			// map "flow number" -> "number of rows in flow"
			Map<Integer, Integer> map = new HashMap<Integer,Integer>();

			checkPortsSize(map, in.ports, inputFlows);
			checkPortsSize(map, wndOut.ports, outputFlows);

			// Composite output including all the rows
			for (int output = 0; output < wndOut.ports.size(); output++)
				if ((group = outputFlows.get(output)) != null)
				{
					List<ErrorRow> indices = errorIndices.get(group);
					Collections.sort(indices);
					for (int i = 0; i < indices.size(); i++)
						indices.get(i).restoreTo(wndOut.ports.get(output));
				}

			// Restore inputs
			for (int input = 0; input < in.ports.size(); input++)
				if ((group = inputFlows.get(input)) != null)
				{
					List<ErrorRow> indices = errorIndices.get(group);
					Collections.sort(indices);
					for (int i = 0; i < indices.size(); i++)
						indices.get(i).restoreTo(in.ports.get(input));
				}

			// Copy attachments from inputs
			for (Integer curOutput : outputFlows.keySet()) 
			{
				for (Integer curInput : inputFlows.keySet()) 
				{
					if (inputFlows.get(curInput).equals(outputFlows.get(curOutput)))
						if (in.ports.size() > curInput)
						{
							DataTable dtIn = in.ports.get(curInput);
							DataTable dtOut = wndOut.ports.get(curOutput);
							if (dtIn.attachment != null && dtOut.attachment == null)
								dtOut.attachment = dtIn.attachment;

							if (dtIn.columnAttachments != null && dtOut.columnAttachments == null)
								dtOut.columnAttachments = dtIn.columnAttachments;
							else
								if (dtIn.columnAttachments != null && dtOut.columnAttachments != null)  // some new values appeared; replacing them with new ones
									for (String column : dtIn.columnAttachments.keySet())
									{
										Map<String, Serializable> ins = dtIn.columnAttachments.get(column);
										Map<String, Serializable> ous = dtOut.columnAttachments.get(column); 
										if(ous != null) // we have something for this column
											for(String key:ins.keySet()) {
												if(!ous.containsKey(key))  // is missed in the output
													ous.put(key, ins.get(key));
											}
										else
											dtOut.columnAttachments.put(column, ins);  //otherwise add it completely
									}

							// Copy also rows' attachments
							// Currently, the attachment is taken from first available input port
							dtIn.reset();
							dtOut.reset();
							while (dtIn.nextRow() && dtOut.nextRow())
								if (dtIn.getCurrentRow().attachments != null) // if there were any attachments in input, copy them to the output
								{
									Long i1 = dtIn.getCurrentRow().attachments == null ? null: (Long)dtIn.getCurrentRow().attachments.get(QSPRConstants.RECORD_ID_ATTACHMENT);
									Long i2 = dtOut.getCurrentRow().attachments == null ? null: (Long)dtOut.getCurrentRow().attachments.get(QSPRConstants.RECORD_ID_ATTACHMENT);
									if( i1 != null && i2 != null && (long)i2 != (long)i1)
										throw new IOException("Data integrity problem - IDS are different: " + i1 + " != " + i2);

									if (dtOut.getCurrentRow().attachments != null)
										dtIn.getCurrentRow().attachments.putAll(dtOut.getCurrentRow().attachments);
									dtOut.getCurrentRow().attachments = dtIn.getCurrentRow().attachments;
								}
						}
				}

			}

			// Logging
			this.out.println(""+in+" -> "+wndOut+"");
			for (Integer curGroup : errorIndices.keySet()) 
			{
				List<ErrorRow> errors = errorIndices.get(curGroup);
				if (errors.size() != 0)
					this.out.println(errors.size()+" errors in input ignored");
			}

			if (currentTask.get().debug > DebugLevel.NONE)
				try{
					Marshaller marshaller = jContext.createMarshaller();
					marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
					marshaller.marshal(wndOut, new File(getDebugDirectory(currentTask.get().id) + "/out-" + currentTask.get().id + "-" + supportedTaskType + ".xml"));
					new ExcelWriter(getDebugDirectory(currentTask.get().id) + "/out-" + currentTask.get().id + "-" + supportedTaskType + ".xls").addWnd(wndOut).finalize();
				}catch(Exception e) {
					e.printStackTrace(System.out);
				}

			this.out.println();
		}catch(Exception e) {
			throw e;
		}
		catch(Throwable ee) {
			throw ee;
		}
		finally {
			Various.molecule = Various.getDefaultCheminfImpl(); // could be required for other purposes in Main Test
		}

		return wndOut;
	}

	WorkflowNodeData calculateInParralel(WorkflowNodeData in, Serializable configuration, int resubmitted) throws Exception
	{
		int suggestedRepostSize = getRepostSize(configuration);

		suggestedRepostSize =  resubmitted == 0? suggestedRepostSize: suggestedRepostSize/(resubmitted + 1);

		if(suggestedRepostSize < 10) suggestedRepostSize = 1;

		out.println("Calculating a task in smaller batches of " + suggestedRepostSize + " rows");
		CalculationTaskSet taskSet = new CalculationTaskSet(this);

		DataTable allMolecules = in.ports.get(0);

		int bags = allMolecules.getRowsSize() / suggestedRepostSize; // number of bags
		int newSize = allMolecules.getRowsSize() / (bags + 1) + 1; // making the bags of a similar size, i.e. avoid to have 30,000 and 1,000 
		suggestedRepostSize = newSize < getRepostSize(configuration) *4/5? newSize : getRepostSize(configuration); // but only if it makes sense!

		in.reset();
		for (int moleculeNumber = 0; moleculeNumber < allMolecules.getRowsSize(); moleculeNumber += suggestedRepostSize)
		{
			WorkflowNodeData wndSmallerBatch = new WorkflowNodeData();

			for (int p = 0; p < in.ports.size(); p++)
			{
				DataTable dt = new DataTable();
				dt.columnAttachments = allMolecules.columnAttachments;
				dt.id = allMolecules.id;
				for (String column : allMolecules.getColumns())
					dt.addColumn(column);
				wndSmallerBatch.addPort(dt);
			}

			for (int currentRecordInFile = 1; currentRecordInFile <= suggestedRepostSize && in.nextRow(); currentRecordInFile++)
			{
				for (int p = 0; p < in.ports.size(); p++)
					wndSmallerBatch.ports.get(p).addRow(in.ports.get(p).getCurrentRow());
			}

			CalculationTask cTask = new CalculationTask();
			cTask.wndInput = wndSmallerBatch;
			cTask.taskName = supportedTaskType;
			cTask.configuration = configuration;
			cTask.setParent(this);
			taskSet.tasks.add(cTask);
		}

		taskSet.allowFailedTasks = true;
		taskSet.post();
		taskSet.calculate(true);

		// Combine results
		DataTable dtResults = new DataTable();
		for (int i = 0; i < taskSet.tasks.size(); i++) //Override with compact row format if necessary
		{
			CalculationTask cTask = taskSet.tasks.get(i);
			if (cTask.error == null)
				if (cTask.getWndOutput().ports.get(0).isCompactRowFormat())
				{
					dtResults = new DataTable(true);
					break;
				}
		}

		for (int i = 0; i < taskSet.tasks.size(); i++)
		{
			CalculationTask cTask = taskSet.tasks.get(i);
			setStatus("Combining results with table no. " + (i + 1) + "...");
			if (cTask.error == null)
				dtResults.addRowsFrom(cTask.getWndOutput().ports.get(0));
			else
				for (int k = 0; k < cTask.size; k++)
				{
					dtResults.addRow();
					dtResults.getCurrentRow().setError(cTask.error);
				}
			taskSet.tasks.set(i, null); // releasing memory for already processed tasks
		}

		setStatus("compacting resulting table and sending results");

		return new WorkflowNodeData(dtResults);
	}

	protected void checkPortsSize(Map<Integer, Integer> numRows, 
			List<DataTable> ports, 
			HashMap<Integer, Integer> mappings) throws Exception
	{
		Integer group;
		Integer size;
		for (int input = 0; input < ports.size(); input++)
			if ((group = mappings.get(input)) != null && ports.get(input) != null)
			{
				if ((size = numRows.get(group)) != null)
				{
					if (size != ports.get(input).getRowsSize())
						throw new Exception(this.supportedTaskType+":Flow datatable size constrain error on port number "+input+" ("+ports.get(input).getRowsSize()+"<>"+size+")");
				} else
					numRows.put(group, ports.get(input).getRowsSize());
			}
	}
	/**
	 * Function control data flows
	 * Indicates which number input ports are under consistency control:
	 *  namely
	 *     1) same number of rows
	 *     2) work with errors -- all ports in the same group will get the same control of errors 
	 *     (i.e., if there is error in one port in a particular row it the same row within the group will be excluded)
	 *     thus the output ports (if they within the same group) will got the same rows marked as error 
	 * @param inputNumber  -- indicates port (i.e., first port is usually data, second port can be e.g. experimental values), starts with 0 
	 * @param groupNumber --  group number (so far not used)
	 */
	protected void setInputFlowGroup(int inputNumber, int groupNumber)
	{
		inputFlows.put(inputNumber, groupNumber);
	}

	protected void setInputFlowGroup(int inputNumber)
	{
		int groupNumber=1;
		inputFlows.put(inputNumber, groupNumber);
	}

	protected void setOutputFlowGroup(int inputNumber, int groupNumber)
	{
		outputFlows.put(inputNumber, groupNumber);
	}


	protected void setOutputFlowGroup(int inputNumber)
	{
		int groupNumber=1;
		outputFlows.put(inputNumber, groupNumber);
	}

	@Override
	public void setParam(String name, String value)
	{
		super.setParam(name, value);
		if ("PIPEIN".equals(name))
			for (int i = 0; i < value.length()/2; i++)
				setInputFlowGroup(Integer.parseInt(""+value.charAt(2*i)), Integer.parseInt(""+value.charAt(2*i+1)));
		if ("PIPEOUT".equals(name))
			for (int i = 0; i < value.length()/2; i++)
				setOutputFlowGroup(Integer.parseInt(""+value.charAt(2*i)), Integer.parseInt(""+value.charAt(2*i+1)));
		params.put(name, value);
	}

	List<Integer> removeStubs(DataTable dt) throws Exception
	{
		List<Integer> indices = new ArrayList<Integer>();

		dt.reset();
		while (dt.nextRow())
			if (WorkflowNodeServer.STUBROW.equals(dt.getCurrentRow().status))
				indices.add(dt.currentRow);
		for (int i = indices.size() - 1; i >= 0; i--)
			dt.deleteRow(indices.get(i));

		return indices;
	}

	void restoreStubs(DataTable dt, List<Integer> indices)
	{
		for (int i = 0; i < indices.size(); i++)
			dt.insertRow(indices.get(i), dt.getStubRow());
	}

	protected long getCalculationTime()
	{
		return (Calendar.getInstance().getTimeInMillis() - taskStartedTime) / 1000;
	}

	/**
	 * Main method to be implemented by all servers
	 * @param task
	 * @param configuration
	 * @return
	 * @throws Exception
	 */

	public abstract WorkflowNodeData calculate(WorkflowNodeData task, Serializable configuration) throws Exception;

	/**
	 * 
	 * Analyze message when calculation failed
	 * If it is critical, e.g. something wrong in data, ctask will not be sent for calculations again
	 * @param message
	 * @return
	 */

	public boolean isCritical(String message) {
		return false;
	}

}

@SuppressWarnings("rawtypes")
class ErrorRow implements Comparable
{
	public int index;
	public String error="";
	public Map<String, Serializable> attachments;

	public ErrorRow(int index)
	{
		this.index = index;
		error = "";
	}

	public ErrorRow addError(String error)
	{
		if(this.error==null)this.error="";
		if(error==null)return this;
		if (!this.error.startsWith(error))
			this.error += error + " ";
		return this;
	}

	@Override
	public boolean equals(Object a)
	{
		return index == ((ErrorRow) a).index;
	}

	public int compareTo(Object a)
	{
		ErrorRow comp = (ErrorRow) a;
		if (comp.index == index)
			return 0;
		else
			return (index > comp.index) ? 1 : -1;
	}

	public ErrorRow setOriginalRow(AbstractDataRow row)
	{
		if (attachments == null)
			attachments = row.attachments;
		if (row.isError())
			addError(row.detailedStatus);
		return this;
	}

	public void restoreTo(DataTable dt)
	{
		AbstractDataRow row = dt.createRow();
		row.setError(error);
		row.attachments = attachments; // FIXME Midnighter: Here we got an error.. we can have wrong attachment here, from which input they are??
		dt.insertRow(index, row);
	}




}