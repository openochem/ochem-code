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

package qspr.tests;

import java.io.IOException;
import java.io.Serializable;

import qspr.metaserver.configurations.ASNNConfiguration;
import qspr.metaserver.configurations.SelectionConfiguration;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.MockTransport;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

/**
 * A mock transport to simulate the model creation process
 * 
 * @author robert
 *
 */
public class ModelCreationMockTransport extends MockTransport
{
	private static final long serialVersionUID = 1L;
	public static double defaultDMValue = 1.0;

	@Override
	public Serializable getTaskResult(Task inTask) throws IOException, ClassNotFoundException
	{
		if (inTask == null)
			inTask = this.inTask;

		WorkflowNodeData wndInput = (WorkflowNodeData) inTask.getData();
		WorkflowNodeData defaultResult = new WorkflowNodeData();
		DataTable dtPredictions = new DataTable(true);
		dtPredictions.addColumn(QSPRConstants.PREDICTION_RESULT_COLUMN);
		dtPredictions.addColumn(QSPRConstants.DM + ":TESTDM");
		for (int i = 0; i < wndInput.ports.get(0).getRowsSize(); i++)
		{
			dtPredictions.addRow();
			dtPredictions.setValue(0, wndInput.ports.get(1).getValue(i, 0));
			dtPredictions.setValue(1, defaultDMValue);
			dtPredictions.getCurrentRow().attachments = wndInput.ports.get(0).getRow(i).attachments;
		}
		defaultResult.addPort(dtPredictions);
		defaultResult.addPort(new DataTable (new ASNNConfiguration()));
		defaultResult.addPort(dtPredictions.getCopy().setId("descriptors"));
		defaultResult.addPort(new DataTable(new SelectionConfiguration()).setId("selection-configuration"));

		return defaultResult;
	}
};
