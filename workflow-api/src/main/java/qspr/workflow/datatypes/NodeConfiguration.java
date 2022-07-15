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

import java.io.Serializable;

import javax.xml.bind.annotation.XmlAnyElement;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElementRef;
import javax.xml.bind.annotation.XmlElementRefs;

public class NodeConfiguration implements Serializable
{

	private static final long serialVersionUID = 10L;

	@XmlAttribute
	// ID of the node in workflow
	public String id;
	
	@XmlAttribute
	// Allows to override task type of the node
	public String taskType; 
	
	@XmlAttribute
	public boolean skipNode = false; 
	
	@XmlAnyElement(lax = true)
	@XmlElementRefs({
	    @XmlElementRef(name="double-list", type=DoubleList.class),
	     @XmlElementRef(name="string-list", type=StringList.class),
	     @XmlElementRef(name="datatable", type=DataTable.class)
	   })
	public Object configuration;
	
	public NodeConfiguration(String id, Serializable configuration)
	{
		this.id = id;
		this.configuration = configuration;
	}
	
	public NodeConfiguration setTaskType(String taskType)
	{
		this.taskType = taskType;
		return this;
	}
	
	public NodeConfiguration setSkipNode(boolean skipNode)
	{
		this.skipNode = skipNode;
		return this;
	}
	
	public NodeConfiguration()
	{
	}
}
