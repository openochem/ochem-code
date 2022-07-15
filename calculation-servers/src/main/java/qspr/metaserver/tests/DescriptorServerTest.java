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

package qspr.metaserver.tests;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import qspr.metaserver.CalculationServer;
import qspr.metaserver.ServerPool;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.cs.DescriptorsAbstractServer;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.TransportFactory;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.descriptorcache.DescriptorsRepository;
import com.eadmet.descriptorcache.DescriptorsRepositoryFactory;

public class DescriptorServerTest 
{

	@Test(timeout = 180000)
	@TaskTest("Descriptors")
	public void descriptorsTest() throws Exception
	{
		String descType = "TestDesc" + Math.random();
		
		int saveTimeout = TransportFactory.waitingTimeBetweenRequests;
		try
		{
			TransportFactory.waitingTimeBetweenRequests = 500;
			DescriptorsConfiguration conf = new DescriptorsConfiguration();
			conf.addDescriptorType(descType, null, null, true, false);
			
			conf.types.get(0).configuration = null;
			
			DataTable dtMols = new DataTable(false);
			dtMols.addColumn(QSPRConstants.SDF_COLUMN);
			dtMols.addRow();
			dtMols.setValue("1");
			dtMols.getCurrentRow().addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, 1);
			dtMols.addRow();
			dtMols.setValue("2");
			dtMols.getCurrentRow().addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, 2);
			dtMols.addRow();
			dtMols.setValue("-3");
			dtMols.getCurrentRow().addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, 3);
			
			dtMols.addRow();
			dtMols.getCurrentRow().addAttachment(QSPRConstants.EXTERNAL_ID, "EXTERNAL-ID-SAMPLE");
			
			Task task = new Task("Descriptors", conf, new WorkflowNodeData(dtMols));
			task.setUser("Test");
		
			TestDescServ testDescSrv = new TestDescServ();
			testDescSrv.supportedTaskType = descType;
			ServerPool.getInstance().servers.add(testDescSrv);
			CalculationServer server = ServerPool.getInstance().getFreeServer("Descriptors");
			
			server.calculateWrapper(task);

			// Do the initial simulated calculation of three molecules. Two normal mols and one error
			DataTable dtDesc = WorkflowNodeData.fromTask(task).ports.get(0);
			assertEquals(1.0, dtDesc.getValue(0, 0));
			assertEquals(2.0, dtDesc.getValue(1, 1));
			assertTrue(dtDesc.getRow(2).isError());
			assertEquals(3, testDescSrv.calculatedMols);
			assertTrue(dtDesc.getRow(3).isError());
			assertTrue(dtDesc.getRow(3).detailedStatus.contains("No stored"));
			
			// Two normal cache entries, one failure to be re-attempted
			testDescSrv.calculatedMols = 0;
			server.calculateWrapper(task);
			assertEquals(1, testDescSrv.calculatedMols);
			
			// Two normal cache entries, one failure to be re-attempted
			testDescSrv.calculatedMols = 0;
			server.calculateWrapper(task);
			assertEquals(1, testDescSrv.calculatedMols);
			
			// No more re-attempts. The failure is permanent
			testDescSrv.calculatedMols = 0;
			server.calculateWrapper(task);
			assertEquals(0, testDescSrv.calculatedMols);
			
			// One new molecule with another descriptor
			testDescSrv.calculatedMols = 0;
			dtMols.setValue(1, 0, "3");
			task = new Task("Descriptors", conf, new WorkflowNodeData(dtMols)).setUser("Test");
			server.calculateWrapper(task);
			assertEquals(1, testDescSrv.calculatedMols);
			dtDesc = WorkflowNodeData.fromTask(task).ports.get(0);
			assertEquals(3.0, dtDesc.getValue(1, "TEST_DESC0"));
		}
		finally
		{
			ServerPool.getInstance().disableTaskType(descType);
			TransportFactory.waitingTimeBetweenRequests = saveTimeout;
			DescriptorsRepository dRepository = DescriptorsRepositoryFactory.getRepository();
			dRepository.clearCache(descType);
		}
	}
}

class TestDescServ extends DescriptorsAbstractServer
{
	public int calculatedMols = 0;
	public TestDescServ()
	{
		supportedTaskType = "TestDesc";
	}

	@Override
	public WorkflowNodeData calculateDescriptors(WorkflowNodeData task, DescriptorsAbstractConfiguration configuration) throws Exception
	{
		DataTable dtDesc = new DataTable(true);
		DataTable dtMols = task.ports.get(0);
		calculatedMols = dtMols.getRowsSize();
		
		dtMols.reset();
		while (dtMols.nextRow())
		{
			Float val = new Float("" + dtMols.getValue());
			dtDesc.addRow();
			if (val > 0)
				dtDesc.setValue("TEST_DESC" + (Math.round(val) % 3), val);
			else
				dtDesc.getCurrentRow().setError("A test simulated error");
			
		}
		
		return new WorkflowNodeData(dtDesc);
	}

}
