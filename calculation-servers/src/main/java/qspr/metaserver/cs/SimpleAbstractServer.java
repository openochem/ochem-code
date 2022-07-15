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

import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

abstract public class SimpleAbstractServer extends WorkflowNodeServer {

	@Override
	public WorkflowNodeData calculate(WorkflowNodeData wndInput,
			Serializable configuration) throws Exception 
	{
	DataTable dtInputData = wndInput.ports.get(0);
	DataTable dtResults = calculateGeneralTask(dtInputData, configuration);
	return new WorkflowNodeData(dtResults);
	}

	abstract DataTable calculateGeneralTask(DataTable dtInputData,
			Serializable configuration) throws Exception;

}
