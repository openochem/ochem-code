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

import javax.xml.bind.annotation.XmlRootElement;

import qspr.metaserver.protocol.Mergable;
import qspr.workflow.structure.WorkflowStructure;

@XmlRootElement
public class WorkflowConfiguration implements Serializable, Mergable
{
	public static final long serialVersionUID = 100L;

	public static final long SAVE_CONFIG = -1;
	public static final long LOAD_CONFIG = -2;

	public String taskId;
	public Long configurationId;
	public NodesConfiguration nodesConfiguration;
	public WorkflowStructure workflowStructure; 
	public int repostSize = 0; // Allows to configure parallelization of calculations 

	public WorkflowConfiguration(String taskId)
	{
		this.taskId = taskId;
	}

	public WorkflowConfiguration(String taskId, NodesConfiguration nodesConfiguration)
	{
		this.taskId = taskId;
		this.nodesConfiguration = nodesConfiguration;
	}

	public WorkflowConfiguration setRepostSize(int repostSize)
	{
		this.repostSize = repostSize;
		return this;
	}

	public void mergeWith(Mergable override)
	{
		WorkflowConfiguration wcOverride = (WorkflowConfiguration) override;
		if (wcOverride.configurationId != null)
			configurationId = wcOverride.configurationId;
		if (wcOverride.nodesConfiguration != null)
			nodesConfiguration = wcOverride.nodesConfiguration;
	}

	public boolean equals(Object obj)
	{
		if(this == obj)
			return true;
		if((obj == null) || (obj.getClass() != this.getClass()))
			return false;
		WorkflowConfiguration test = (WorkflowConfiguration)obj;

		return (equals(configurationId, test.configurationId) && equals(taskId, test.taskId));
	}

	public boolean equals (Object a, Object b)
	{
		if (a == null)
			return (b == null);
		return a.equals(b);
	}

	public WorkflowConfiguration setNodesConfiguration(NodesConfiguration configuration)
	{
		nodesConfiguration = configuration;
		return this;
	}

	public WorkflowConfiguration addNodeConfiguration(NodeConfiguration configuration)
	{
		if (nodesConfiguration == null)
			nodesConfiguration = new NodesConfiguration();
		nodesConfiguration.nodes.add(configuration);

		return this;
	}

	public WorkflowConfiguration()
	{

	}

}
