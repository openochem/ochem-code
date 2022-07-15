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
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.net.MalformedURLException;
import java.util.List;
import java.util.concurrent.TimeoutException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Assert;
import org.junit.Test;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.entities.Molecule;
import qspr.entities.ReadyModelAttachment;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.configurations.BalloonConfiguration;
import qspr.metaserver.configurations.ModelAbstractDataDrivenConfiguration;
import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.CSTransport;
import qspr.metaserver.transport.Transport;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.configurations.CDSModelData;
import qspr.util.MoleculePeer;
import qspr.util.NCBI_Utility;
import qspr.util.WrapperThread;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;

@PeriodicTest(schedule = ScheduleType.HOURLY, type = "general")
public class PeriodicTests
{
	private static transient final Logger logger = LogManager.getLogger(PeriodicTests.class);

	public String[] testSDFs;

	public static final String ETHYLACETATE = "8857\n  -OEChem-12090808052D\n\n 14 13  0     0  0  0  0  0  0999 V2000\n    3.7321    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5981    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.4641    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.1996   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.9966   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.7741    0.2131    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1541    1.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3100    1.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4631    1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6900    0.2131    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  3  1  0  0  0  0\n  1  5  1  0  0  0  0\n  2  5  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  7  1  0  0  0  0\n  3  8  1  0  0  0  0\n  4  9  1  0  0  0  0\n  4 10  1  0  0  0  0\n  4 11  1  0  0  0  0\n  5  6  1  0  0  0  0\n  6 12  1  0  0  0  0\n  6 13  1  0  0  0  0\n  6 14  1  0  0  0  0\nM  END";
	public static final String ALANINE = "Alanine\n  Marvin  12100811572D\n\n  6  5  0  0  0  0            999 V2000\n   -0.6811    0.9144    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0333    1.3269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7478    0.9144    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0333    2.1519    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.6811    0.0894    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.3956    1.3269    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n  2  4  1  0  0  0  0\n  1  5  1  0  0  0  0\n  1  6  1  0  0  0  0\nM  END";
	public static final String BENZOCAINE = "Benzocaine\n  Marvin  12100812302D\n\n 12 12  0  0  0  0            999 V2000\n   -0.2501    2.8409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.2501    2.0159    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9645    1.6034    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9645    0.7784    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6790    2.0159    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6790    2.8409    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9645    3.2534    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.9645    4.0784    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -1.6790    4.4909    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.2501    4.4909    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.2501    5.3159    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.4644    5.7284    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  5  6  2  0  0  0  0\n  6  7  1  0  0  0  0\n  7  8  1  0  0  0  0\n  1  7  2  0  0  0  0\n  8  9  2  0  0  0  0\n  8 10  1  0  0  0  0\n 10 11  1  0  0  0  0\n 11 12  1  0  0  0  0\nM  END";


	public PeriodicTests() throws IOException, TimeoutException
	{
		String[] smiles = new String[]{"CCC", "CCCC", "CCCCC", "CCCCCC"};
		testSDFs = new String[smiles.length];
		for (int i = 0; i < smiles.length; i++)
			testSDFs[i] = Various.molecule.convertToSDFUnblocked(smiles[i]).replace("$$$$", "").trim();
	}

	@Test(timeout = 20000)
	@NotifyOnFailure(minConsequentFailures = 3)
	public void metaserverTest() throws Exception
	{
		Transport transport = new CSTransport();
		transport.executeCommand(new Command(Command.CL_QUERY_TASK, 666).sid("Test runner"));
	}

	@Test(timeout = 120000)
	@NotifyOnFailure(minConsequentFailures = 3)
	public void conversion3DTest() throws Exception
	{
		checkIfTaskSupported(QSPRConstants.BALLOON);
		StructureOptimisationConfiguration conf = new BalloonConfiguration();
		conf.bypassCache = true;

		Task task = new Task(conf.getTaskType(), conf, new WorkflowNodeData(new DataTable(Various.molecule.convertToSDFUnblocked("CC(=O)"))), false);
		task = calculateTask(task);
	}

	/**
	 * A simple test to check whether the Xemistry indexing tables are OK
	 */
	@Test(timeout = 10000)
	@NotifyOnFailure(minConsequentFailures = 2)
	public void xemistryTablesStatus() throws Exception
	{
		WrapperThread t = new WrapperThread()
		{
			@Override
			public void wrapped() throws Exception
			{
				Globals.session().createSQLQuery("select * from StructureQuery limit 1").list();
			}
		};

		t.run();

		if (t.exception != null)
			throw t.exception;
	}

	private Task calculateTask(Task inTask) throws Exception
	{
		inTask.setUser(QSPRConstants.TEST_USER_PREFIX+"User");
		Task task = new CalculationClient().setSid("Test runner").calculateTask(inTask);
		task.check();
		return task;
	}

	@Test(timeout = 120000)
	@NotifyOnFailure(minConsequentFailures = 3)
	public void amesModelTest() throws Exception
	{
		WrapperThread wrapper = new DataDrivenTestWrapper("Ames-Test") {

			@Override
			public void wrappedTest() throws Exception 
			{
				ModelApplier applier = new ModelApplier();
				applier.compoundsProvider.basket.addMolecule(MoleculePeer.getMolecule(ALANINE));
				applier.compoundsProvider.basket.addMolecule(MoleculePeer.getMolecule(ETHYLACETATE));

				applier.addModel(Repository.model.getByPublicId(1));
				applier.useCache = false;

				applier.startAndWait();

				Assert.assertTrue(!applier.isError());
				Assert.assertTrue(!applier.modelTasks.get(0).isError());
				DataTable dtPredictions = applier.modelTasks.get(0).wndResult.ports.get(0);
				Assert.assertEquals(0, dtPredictions.errorCount());
				Assert.assertEquals(2, dtPredictions.getRowsSize());
			}
		};
		wrapper.run();

		if (wrapper.exception != null)
			throw wrapper.exception;
	}





	@Test(timeout = 10000)
	@NotifyOnFailure(minConsequentFailures = 3)
	public void used_memoryTest() throws Exception
	{
		final int MAXMEMORY = 5000;
		try
		{
			MemoryPoolMXBean pool = getPool();
			log("Memory test for pool "+pool.getName()+" of type "+pool.getType());

			long maxMemory = Math.round(1.0 * maxUsedMemory() / 1024 / 1024); // In MB
			long usedMemory = Math.round(1.0 * usedMemory() / 1024 / 1024);

			log("Max used memory: " + maxMemory + " Mb");
			log("Currently used memory:" + usedMemory + " Mb");

			if (usedMemory > MAXMEMORY)
				throw new Exception("Currently used memory is " + usedMemory + "MB, Limit is " + MAXMEMORY + "MB" );

			if (maxMemory > MAXMEMORY)
				throw new Exception("Maximal used memory since last test run was " + maxMemory + "MB, Limit is " + MAXMEMORY + "MB");
		}
		finally
		{
			resetMaxUsedMemory();
		}
	}

	@Test(timeout = 120000)
	public void ncbiSearchTest() throws Exception
	{
		Globals.startAllTransactions();
		List<Molecule> mols = NCBI_Utility.getMoleculeByName("viagra", 1);	
		Globals.rollbackAllTransactions();
		Assert.assertEquals(1, mols.size());
	}


	protected static MemoryPoolMXBean getPool()
	{
		List<MemoryPoolMXBean> pools = ManagementFactory.getMemoryPoolMXBeans();
		MemoryPoolMXBean result = null;
		for(MemoryPoolMXBean pool : pools) 
			if (result == null || (result.getUsage().getMax() < pool.getUsage().getMax()))
				result = pool;
		return result;
	}

	public static long usedMemory()
	{
		MemoryPoolMXBean pool = getPool();
		return pool.getUsage().getUsed();
	}

	public static long maxUsedMemory()
	{
		MemoryPoolMXBean pool = getPool();
		return pool.getPeakUsage().getUsed();
	}

	public static void resetMaxUsedMemory()
	{
		MemoryPoolMXBean pool = getPool();
		pool.resetPeakUsage();
	}

	private void log(String msg)
	{
		logger.info("[Periodic tests] " + msg);
	}

	static public void checkIfTaskSupported(String taskType) throws MalformedURLException, IOException, ClassNotFoundException
	{
		if (!new CalculationClient("PeriodicTests").getSupportedTaskTypes().contains(taskType))
			throw new UserFriendlyException("Task " + taskType + " is not supported by MetaServer");
	}

	public static void main(String[] args) throws MalformedURLException, IOException, ClassNotFoundException {
		new WrapperThread()
		{

			@SuppressWarnings("unchecked")
			@Override
			public void wrapped() throws Exception
			{
				qspr.entities.Attachment<ReadyModelAttachment> a = (qspr.entities.Attachment<ReadyModelAttachment>) Globals.session().get(qspr.entities.Attachment.class, 605885L);
				ReadyModelAttachment rma = a.getObject();
				CDSModelData cds = (CDSModelData) rma.modelData;
				System.out.println(((ModelAbstractDataDrivenConfiguration)cds.methodSpecificData).trainingSet.get());
			}
		}.run();
	}
}
