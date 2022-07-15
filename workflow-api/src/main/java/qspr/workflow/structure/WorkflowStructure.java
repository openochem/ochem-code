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

package qspr.workflow.structure;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;


@XmlRootElement(name = "workflow")
public class WorkflowStructure implements Serializable
{
	private static final long serialVersionUID = 1L;

	@XmlAttribute
	public Integer insockets = 1;

	@XmlAttribute
	public Integer outsockets = 1;

	@XmlElementWrapper(name = "nodes")
	@XmlElement(name = "node")
	public List<StructureNode> nodes;

	@XmlElementWrapper(name = "connections")
	@XmlElement(name = "connection")
	public List<Connection> connections;	
	
	private StructureNode lastNode;
	private Connection lastConnection;
	
	public static WorkflowStructure fromXml(JAXBContext jContext, String filename) throws Exception
	{
		Unmarshaller ml = jContext.createUnmarshaller();
		return (WorkflowStructure)ml.unmarshal(new File(filename));
	}
	
	public void addNode(String nodeId, String taskName)
	{
		StructureNode node = new StructureNode();
		node.id = nodeId;
		node.taskName = taskName;
		
		
		String sourceNode = lastNode != null ? lastNode.id : "in";
		
		if (nodes == null)
			nodes = new ArrayList<StructureNode>();
		
		nodes.add(node);
		
		if (connections == null)
			connections = new ArrayList<Connection>();
		
		connections.add(new Connection(sourceNode, node.id));
		
		if (lastConnection == null)
		{
			lastConnection = new Connection(node.id, "out");
			connections.add(lastConnection);
		} else
			lastConnection.from = node.id;
		
		lastNode = node;
	}
}


