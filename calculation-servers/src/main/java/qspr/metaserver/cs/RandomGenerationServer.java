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

import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

public class RandomGenerationServer extends WorkflowNodeServer
{
	public RandomGenerationServer()
	{
		supportedTaskType = DescriptorsConfiguration.RANDOM;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}
	
	@Override
	public WorkflowNodeData calculate(WorkflowNodeData task,
			Serializable receivedConfiguration) throws Exception 
	{
		DataTable dtResult = new DataTable(true);
		double num = (Double) receivedConfiguration;
		
		
		for (int i = 0; i < num; i++)
			dtResult.addColumn("Random+" + (i + 1));
		
		for (int i = 0; i < task.ports.get(0).getRowsSize(); i++)
		{
			dtResult.addRow();
			for (int k = 0; k < num; k++)
				dtResult.setValue(k, Double.valueOf(Math.random()));
		}
		
		
		return new WorkflowNodeData(dtResult);
	}
}
