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

import java.io.Serializable;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsEmptyConfiguration;
import qspr.workflow.datatypes.WorkflowNodeData;

// Is there really any value of this class?
public abstract class DescriptorsAbstractServer extends WorkflowNodeServer {

	@Override
	public WorkflowNodeData calculate(WorkflowNodeData wndInput,
			Serializable configuration) throws Exception {

		DescriptorsAbstractConfiguration configurationDesc = (DescriptorsAbstractConfiguration) configuration;

		if (configurationDesc == null)
			configurationDesc = new DescriptorsEmptyConfiguration();

		return calculateDescriptors(wndInput, configurationDesc);
	}

	/**
	 * General interface to calculate Descriptor Task
	 * 
	 * @param task
	 * @param configuration
	 * @return
	 * @throws Exception
	 */
	@Override
	protected int getRepostSize(Serializable configuration)
	{
		if (configuration instanceof DescriptorsAbstractConfiguration)
			if (((DescriptorsAbstractConfiguration) configuration).repostSize != null)
				return ((DescriptorsAbstractConfiguration) configuration).repostSize;
		
		return repostSize;
	}
	
	public abstract WorkflowNodeData calculateDescriptors(
			WorkflowNodeData task, DescriptorsAbstractConfiguration configuration)
			throws Exception;

}
