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
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Assert;
import org.junit.Test;

import qspr.exceptions.CalculationException;
import qspr.metaserver.ServerPool;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.LocalTransport;
import qspr.metaserver.transport.TransportFactory;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

public class WorkflowNodeServerTest
{
	final String SERVERNAME = "TestTask";
	
	private static transient final Logger logger = LogManager.getLogger(WorkflowNodeServerTest.class);

	public void parralelCalculationTest(Double[] data, boolean compact) throws Exception
	{
		logger.info("Started " + SERVERNAME + " ...");
		TransportFactory.defaultTransportClass = LocalTransport.class;
		TransportFactory.waitingTimeBetweenRequests = 10;
		LocalTransport.allowQueue = true;
		final boolean comp = compact;
		WorkflowNodeServer server = new WorkflowNodeServer()
		{
			@Override
			public WorkflowNodeData calculate(WorkflowNodeData task, Serializable configuration) throws Exception
			{
				int seed = ((Double)task.ports.get(0).getValue(0, 0)).intValue();
				Random r = new Random(seed);
				if (((Double)task.ports.get(0).getValue(0, 0)) < 0)
					throw new CalculationException("Simulated failure");
				
				DataTable dtRes = new DataTable(comp);
				double rand = r.nextInt(1000);
				
				// One identity column and one random one, for testing columns merging
				dtRes.addColumn(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM);
				dtRes.addColumn("" + rand);
				
				dtRes.addRow();
				dtRes.setValue(0, 0, task.ports.get(0).getValue(0, 0));
				dtRes.setValue(0, 1, rand);
				
				logger.info(dtRes.toStringColumns());
				
				return new WorkflowNodeData(dtRes);
			}
		};
		
		server.supportedTaskType = SERVERNAME;
		server.repostSize = 1; // the server calculates each input sample in a separate run
		
		ServerPool.getInstance().servers.clear();
		ServerPool.getInstance().servers.add(server);
		
		Task task = new Task(SERVERNAME, "Configuration", new WorkflowNodeData(new DataTable(data, "Input")));
		server.calculate(task);
		
		task.check();
		DataTable dtOut = WorkflowNodeData.fromTask(task).ports.get(0);
		
		// Check if the results have been properly merged
		for (int i=0; i<data.length; i++)
			if (data[i] < 0)
				Assert.assertTrue(dtOut.getRow(i).isError());
			else
				Assert.assertEquals(data[i], dtOut.getValue(i, 0));
		
		// Check multi-columns merge
		Set<String> ucols = new HashSet<String>();
		for (int i=0; i<data.length; i++)
			if (data[i] > 0)
			{
				ucols.add(""+data[i]);
				Assert.assertEquals(dtOut.getValue(i, ucols.size()).toString(), dtOut.getColumn(ucols.size()));
			}
		
		Assert.assertEquals(ucols.size() + 1, dtOut.getColumnsSize());
	}
	
//	@Test(timeout = 20000)
	public void parralelCalculationTestDifferentColumnsCompact() throws Exception
	{
		parralelCalculationTest(new Double[]{1.0, 2.0, -3.0, 4.0}, true);
	}
	
//	@Test(timeout = 20000)
	public void parralelCalculationTestSameColumnsCompact() throws Exception
	{
		parralelCalculationTest(new Double[]{1.0, 1.0, -1.0, 1.0}, true);
	}
	
//	@Test(timeout = 20000)
	public void parralelCalculationTestDifferentColumns() throws Exception
	{
		parralelCalculationTest(new Double[]{1.0, 2.0, -3.0, 4.0}, false);
	}
	
//	@Test(timeout = 20000)
	public void parralelCalculationTestSameColumns() throws Exception
	{
		parralelCalculationTest(new Double[]{1.0, 1.0, -1.0, 1.0}, false);
	}

	@Test(timeout = 20000)
	public void noTest() throws Exception
	{
		Assert.assertEquals( 1, 1);

	}

}