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

import java.io.File;
import java.io.InputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;

import qspr.configuration.JAXBContextFactory;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.cs.WorkflowNodeServer;
import qspr.metaserver.protocol.Task.DebugLevel;
import qspr.workflow.datatypes.NodeConfiguration;
import qspr.workflow.datatypes.NodesConfiguration;
import qspr.workflow.datatypes.WorkflowConfiguration;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.structure.Connection;
import qspr.workflow.structure.StructureNode;
import qspr.workflow.structure.WorkflowStructure;
import qspr.workflow.utils.ExcelWriter;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.ResourceLoader;

public class Workflow extends WorkflowNodeServer 
{
	static public final JAXBContext jContext;

	//private String workflowsDirectory = null;

	static 
	{
		try {
			jContext = JAXBContextFactory.get("qspr.workflow.datatypes:qspr.workflow.structure:qspr.metaserver.configurations");
		} catch (JAXBException e) {
			e.printStackTrace();
			throw new ExceptionInInitializerError(e);
		}
	}

	public MyCalculationClient myClient = new MyCalculationClient();
	public List<Node> nodes = new ArrayList<Node>();
	Map<String, Node> nodesMap = new HashMap<String, Node>();

	public SocketSet inputSet;
	public SocketSet outputSet;

	public Workflow()
	{
		supportedTaskType = QSPRConstants.Workflow;
	}

	public void setParam(String name, String value)
	{
		super.setParam(name, value);
	}

	@Override
	protected int getRepostSize(Serializable configuration)
	{
		if (configuration instanceof WorkflowConfiguration)
			return ((WorkflowConfiguration) configuration).repostSize;
		return repostSize;
	}


	private void loadStructure(WorkflowStructure structure) throws Exception
	{

		inputSet = new SocketSet(structure.insockets, true);
		outputSet = new SocketSet(structure.outsockets, false);

		nodes.clear();
		nodesMap.clear();

		for (StructureNode snode : structure.nodes) 
		{
			if (nodesMap.get(snode.id) != null)
				throw new Exception("Error in workflow descriprion - duplicate nodes with ID '"+snode.id+"'");

			Node node = new Node();
			node.id = snode.id;
			node.taskName = snode.taskName;
			//node.configuration = snode.configuration;
			node.input = new SocketSet(snode.insockets, false);
			node.output = new SocketSet(snode.outsockets, true);
			node.parent = this;
			nodes.add(node);
			nodesMap.put(node.id, node);
		}

		for (Connection con : structure.connections) 
		{			
			SocketSet socketsTo = null;
			SocketSet socketsFrom = null;

			if (con.from.equals("in"))
				socketsFrom = inputSet;
			else
				socketsFrom = nodesMap.get(con.from).output;

			if (con.to.equals("out"))
				socketsTo = outputSet;
			else
				socketsTo = nodesMap.get(con.to).input;

			socketsTo.sockets.set(con.toPort, socketsFrom.sockets.get(con.fromPort));
		}

		for (Node node : nodes)
			node.input.checkInputs();

		return;
	}

	private void loadConfiguration(WorkflowConfiguration wc, WorkflowNodeData task) throws Exception
	{
		String workflowsDirectory = "classpath:com/eadmet/workflows";
		WorkflowStructure structure = null;

		if (wc.workflowStructure == null)
		{
			String path =  workflowsDirectory + "/" +wc.taskId.toLowerCase()+".xml";
			InputStream is = ResourceLoader.getInputStreamByPath(path);
			if (is != null)
			{
				structure = (WorkflowStructure) jContext.createUnmarshaller().unmarshal(is);
				out.println("Loading structure from "+path);
				is.close();
			}
		}
		else
			structure = wc.workflowStructure;

		if (structure == null)
			throw new Exception("There is no such workflow: " + wc.taskId);

		loadStructure(structure);


		if (wc.nodesConfiguration != null)
			loadNodesConfiguration(wc.nodesConfiguration);
	}

	private void loadNodesConfiguration(NodesConfiguration conf)
	{
		for (NodeConfiguration node : conf.nodes) 
		{		
			Node realNode = this.nodesMap.get(node.id);
			if (realNode != null)
			{
				realNode.configuration = node.configuration;
				realNode.skip = node.skipNode;
				// We can override the task type
				if (node.taskType != null)
					realNode.taskName = node.taskType;
			}
		}		
	}

	private WorkflowNodeData processSystemTask(Long taskCode, WorkflowNodeData task) throws Exception
	{
		throw new Exception("Unsupported system task");
	}

	@Override
	public WorkflowNodeData calculate(WorkflowNodeData task, Serializable configuration) throws Exception
	{
		myClient.sid = QSPRConstants.Workflow;
		myClient.out = out;
		try
		{	
			WorkflowConfiguration wc = (WorkflowConfiguration) configuration;
			if (wc.taskId.equals("SYSTEM"))
				return processSystemTask(wc.configurationId, task);

			loadConfiguration(wc, task);
			wc = null;
			currentTask.get().clearConfig();

			if (task.ports.size() != inputSet.sockets.size())
				out.println("WARNING: Number of inputs provided ("+task.ports.size()+") is different from workflow's input size (" + inputSet.sockets.size() + ")");

			for (int i = 0; i < Math.min(task.ports.size(), inputSet.sockets.size()); i++)
				inputSet.setValue(i, task.ports.get(i));

			// Fill not provided inputs with empty
			for (int i = 0; i < inputSet.sockets.size(); i++)
				if (i >= task.ports.size() || task.ports.get(i) == null)
				{
					// Allow nulls for workflow input
					out.println("Using null datatable for input " + i);
					inputSet.setValue(i, new NullDataTable());
				}

			// Do not fail if the metaserver is currently down
			myClient.setTolerateMetaserverDown();


			boolean nochEinMalBitte = true;
			while (nochEinMalBitte)
			{
				nochEinMalBitte = false;
				for (int nodeNum = 0; nodeNum < nodes.size(); nodeNum++)
				{
					Node node = nodes.get(nodeNum);
					if (!node.ready && node.readyForCalculation())
					{

						Marshaller marshaller = Workflow.jContext.createMarshaller();
						marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, new Boolean(true));

						nochEinMalBitte = true;
						setStatus(myClient.originalStatus = "Processing task "+node.taskName);
						if (currentTask.get().debug > DebugLevel.NONE)
						{
							marshaller.marshal(node.input.getWorkflowNodeData(), new File(getDebugDirectory(currentTask.get().id) + "/in-" + currentTask.get().id + "-" + node.id + ".xml"));
							new ExcelWriter(getDebugDirectory(currentTask.get().id) + "/in-" + currentTask.get().id + "-" + node.id + ".xls").addWnd(node.input.getWorkflowNodeData()).finalize();
						}

						node.calculate();

						if (currentTask.get().debug > DebugLevel.NONE)
						{
							marshaller.marshal(node.output.getWorkflowNodeData(), new File(getDebugDirectory(currentTask.get().id) + "/out-" + currentTask.get().id + "-" + node.id + ".xml"));
							new ExcelWriter(getDebugDirectory(currentTask.get().id) + "/out-" + currentTask.get().id + "-" + node.id + ".xls").addWnd(node.output.getWorkflowNodeData()).finalize();
						}

					}
				}
			}
			return outputSet.getWorkflowNodeData();

		} catch (Exception e)
		{
			e.printStackTrace();
			throw e;
		}	
	}

	class MyCalculationClient extends CalculationClient
	{
		public String originalStatus = null;

		public void setStatus(String taskStatus)
		{
			super.setStatus(taskStatus);
			String newStatus = originalStatus;
			if (taskStatus != null && !"".equals(taskStatus))
				newStatus += " - "+taskStatus;
			if (!newStatus.equals(Workflow.this.getStatus()))
				Workflow.this.setStatus(newStatus);
		}
	}

}
