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
import static org.junit.Assert.fail;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Map;
import java.util.Random;

import javax.xml.bind.JAXBContext;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;

import qspr.OCHEMConfiguration;
import qspr.configuration.JAXBContextFactory;
import qspr.dao.ChemInfEngine;
import qspr.dao.Various;
import qspr.metaserver.CalculationServer;
import qspr.metaserver.ServerPool;
import qspr.metaserver.configurations.*;
import qspr.metaserver.configurations.DeepChemConfiguration.DeepChemMethod;
import qspr.metaserver.configurations.DescriptorsMORDREDConfiguration.Descr;
import qspr.metaserver.configurations.DescriptorsSIRMSConfiguration.Labeling;
import qspr.metaserver.configurations.KGCNNConfiguration.KGCNN;
import qspr.metaserver.configurations.SKLConfiguration.SKLMethod;
import qspr.metaserver.configurations.LabelWeighting.ClassWeighting;
import qspr.metaserver.configurations.LabelWeighting.PropertyWeighting;
import qspr.metaserver.configurations.ModelAbstractConfiguration.ScalingType;
import qspr.metaserver.cs.ASNNServer;
import qspr.metaserver.cs.CrossValidationServer;
import qspr.metaserver.cs.WorkflowNodeServer;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.CSTransport;
import qspr.metaserver.transport.LocalTransport;
import qspr.metaserver.transport.Transport;
import qspr.metaserver.transport.TransportFactory;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowConfiguration;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.MemoryUtils;

// This a test, run automatically from MultiServer upon update

public class MainTest {
	private static transient final Logger logger = LogManager.getLogger(MainTest.class);
	private static String PATH = "files/test-samples/datatables/";

	final static int MULTIPLIER = 1;

	static JAXBContext jc;
	static Map<String, Serializable> map;

	private static WorkflowNodeData wnd_molecules = null;
	private static DataTable dt_ALOGPS, dt_PYMOL, dt_mopac, dt_default;

	boolean longDragon = false;

	@BeforeClass
	public static void init() throws Exception {
		jc = JAXBContextFactory.get("qspr.workflow.datatypes:qspr.workflow.structure");
		LocalTransport.allowQueue = true;
		TransportFactory.defaultTransportClass = LocalTransport.class;

		Various.molecule = Various.getDefaultCheminfImpl();

		if(Various.molecule == null) throw new UserFriendlyException("Chemoinformatics engine was not enabled (check version-template.xml), the default is : " + OCHEMConfiguration.getCheminfEngine());

		map = Collections.singletonMap(QSPRConstants.CHEMENGINE, "" + Various.molecule.engine);
		dt_ALOGPS = DataTable.fromXml(jc, PATH + "testALOGPS.xml");
		dt_PYMOL = DataTable.fromXml(jc, PATH + "testPYMOL.xml");
		dt_mopac = DataTable.fromXml(jc, PATH + "dt_test_mopacboinc_molecule.xml");
		dt_default = DataTable.fromXml(jc, PATH + "dt_molecules.xml");

		dt_ALOGPS.setInfEngine(Various.molecule.engine);
		dt_PYMOL.setInfEngine(Various.molecule.engine);
		dt_mopac.setInfEngine(Various.molecule.engine);
		dt_default.setInfEngine(Various.molecule.engine);

		wnd_molecules = new WorkflowNodeData(dt_default);
	}

	@AfterClass
	public static void finaliz() {
		TransportFactory.defaultTransportClass = CSTransport.class;
		LocalTransport.allowQueue = false;
		ServerPool.getInstance().disableTaskType("Increment");
		ServerPool.getInstance().disableTaskType("TestModel");
	}

	@Test(timeout = MULTIPLIER * 180000)
	@TaskTest(QSPRConstants.MolStandartizer)
	public void molStandartisationTest() throws Exception {
		StandartizationOptions options = new StandartizationOptions();
		ChemInfEngine standartizer = Various.defaultEngine;
		options.setDefaults(standartizer);
		runStandardizationTestWithOptions(options);
	}

	public void runStandardizationTestWithOptions(StandartizationOptions options) throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.MolStandartizer);
		Task task;
		server.calculateWrapper(
				task = new Task(QSPRConstants.MolStandartizer, options, new WorkflowNodeData(dt_default.getDeepCopy())));
		DataTable dtResponse = WorkflowNodeData.fromTask(task).ports.get(0);
		assertEquals(0, dtResponse.errorCount());
	}

	@Test(timeout = MULTIPLIER * 50000)
	@TaskTest(QSPRConstants.BAGGING)
	public void baggingTest() throws Exception {
		System.out.println(MemoryUtils.memorySummary());

		BaggingConfiguration bc = new BaggingConfiguration();
		bc.keepIndividualPredictions = true;
		bc.ensembleSize = 6;
		bc.taskConfiguration = new ASNNConfiguration();
		bc.taskName = QSPRConstants.ASNN;
		//bc.mixtureValidation = MixtureValidation.MIXTURE;

		DataTable descriptors = generateDescriptorsDatatable(100, 20);
		DataTable values = generateLabelsDatatable(100);
		WorkflowNodeData wnd = new WorkflowNodeData();
		wnd.addPort(descriptors);
		wnd.addPort(values);

		DataTable results = generateLabelsDatatable(200);
		DataTable model = new DataTable(new ASNNConfiguration());
		WorkflowNodeData rwnd = new WorkflowNodeData();
		rwnd.addPort(results);
		rwnd.addPort(model);

		int oldWaitingTime = TransportFactory.waitingTimeBetweenRequests;
		Transport oldTransport = TransportFactory.getThreadTransport();

		try {
			BaggingMockTransport mt = new BaggingMockTransport(rwnd);
			TransportFactory.setThreadTransport(mt);
			TransportFactory.waitingTimeBetweenRequests = 50;

			DataTable result = runTest(QSPRConstants.BAGGING, bc, wnd);

			assertEquals(result.getRowsSize(), 100);
		} finally {
			TransportFactory.setThreadTransport(oldTransport);
			TransportFactory.waitingTimeBetweenRequests = oldWaitingTime;
		}
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(QSPRConstants.DESCRIPTORS)
	public void trivialDescriptorsTest() throws Exception {
		assertTrue(true);
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(DescriptorsConfiguration.QNPR)
	public void QNPRTest() throws Exception {
		DescriptorsQNPRConfiguration c = new DescriptorsQNPRConfiguration();
		DataTable dtResult = runTest(c, wnd_molecules);
		assertTrue(dtResult.getColumnsSize() > 0);
		assertTrue(dtResult.getRowsNoErrorsSize() == dtResult.getRowsSize());
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(DescriptorsConfiguration.SIRMS)
	public void SIRMSTest() throws Exception {
		DescriptorsSIRMSConfiguration sirms = new DescriptorsSIRMSConfiguration(1, 4);
		sirms.descriptorTypes.add(Labeling.elm);
		sirms.descriptorTypes.add(Labeling.REFRACTIVITY);
		sirms.descriptorTypes.add(Labeling.CHARGE);
		sirms.descriptorTypes.add(Labeling.HB);
		sirms.descriptorTypes.add(Labeling.LOGP);
		DataTable dtResult = runTest(sirms, wnd_molecules);
		assertTrue(dtResult.getColumnsSize() > 0);
	}


	@Test(timeout = MULTIPLIER * 5 * 15000)
	@TaskTest(DescriptorsConfiguration.EPA)
	public void EPATest() throws Exception {
		DescriptorsAbstractConfiguration epa = (DescriptorsAbstractConfiguration) Class.forName("qspr.metaserver.configurations.DescriptorsEPAConfiguration").newInstance();
		DataTable dtResult = runTest(epa, wnd_molecules);
		if(dtResult.errorCount() > 0) {
			String errors = "";
			for(int i =0;i<dtResult.getRowsSize();i++)
				if(dtResult.getRow(i).isError())errors += dtResult.getRow(i).detailedStatus+ "\n";
			fail(errors);
		}
		assertTrue(dtResult.getColumnsSize() > 0);
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(DescriptorsConfiguration.OEstate)
	public void OEstateTest() throws Exception {
		DataTable dtResult = runTest(new DescriptorsOEstateConfiguration(), wnd_molecules);
		assertTrue(dtResult.containsColumn("MW"));
	}


	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(DescriptorsConfiguration.ALogPS)
	public void alogpsTest() throws Exception {
		DataTable dtResult = runTest(DescriptorsConfiguration.ALogPS, new DescriptorsEmptyConfiguration(), wnd_molecules);
		assertTrue(dtResult.containsColumn("ALogPS_logP"));
	}

	@Test(timeout = MULTIPLIER * 2 * 15000)
	@TaskTest(DescriptorsConfiguration.MAP4)
	public void map4Test() throws Exception {
		DataTable dtResult = runTest(new DescriptorsMAP4Configuration(), wnd_molecules);
		assertTrue(dtResult.containsColumn("map4_1"));
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(DescriptorsConfiguration.MOLD2)
	public void mold2Test() throws Exception {
		DataTable dtResult = runTest(new DescriptorsMOLD2Configuration(), wnd_molecules);
		assertTrue(dtResult.containsColumn("D001"));
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(DescriptorsConfiguration.Spectrophores)
	public void SpectrophoresTest() throws Exception {
		DataTable dtResult = runTest(new DescriptorsSpectrophoresConfiguration(), wnd_molecules);
		assertTrue(dtResult.containsColumn("SpectrophoresPartial_1"));
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(QSPRConstants.SELECTION)
	public void selectionTest() throws Exception {
		DataTable dtData = new DataTable(true);
		dtData.addColumns(5);
		for (int sample = 0; sample < 100; sample++) {
			double a = Math.random();
			dtData.addRow();
			dtData.setValue(0, new Double(a));
			dtData.setValue(1, new Double(a + (Math.random() - 0.5) / 20));
			dtData.setValue(2, new Double(a + (Math.random() - 0.5) / 20));
			dtData.setValue(3, new Double(1));
			dtData.setValue(4, Math.random());
		}

		SelectionConfiguration conf = new SelectionConfiguration();
		conf.correlationThreshold = 0.5;
		conf.numDifferentValues = 3;

		DataTable dtResult = runTest(QSPRConstants.SELECTION, conf, dtData);
		assertEquals(2, dtResult.getColumnsSize());
	}

	@Test(timeout = MULTIPLIER * 10000)
	@TaskTest(QSPRConstants.CROSSVALIDATION)
	public void crossValidationTest() throws Exception {
		for (int i = 1; i <= 3; i++)
			ServerPool.getInstance().servers.add(new TestModelServ());
		int saveInterval = TransportFactory.waitingTimeBetweenRequests;
		try {
			TransportFactory.waitingTimeBetweenRequests = 50;
			CrossValidationServer server = (CrossValidationServer) ServerPool.getInstance()
					.getFreeServer(QSPRConstants.CROSSVALIDATION);
			server.allowLocalCalculations = true;
			CrossValidationConfiguration cvConf = new CrossValidationConfiguration();
			cvConf.ensembleSize = 2;
			cvConf.taskName = "TestModel";
			cvConf.taskConfiguration = new ASNNConfiguration();
			((ModelAbstractConfiguration) cvConf.taskConfiguration).saveModels = false;
			//cvConf.mixtureValidation = MixtureValidation.MIXTURE;
			Task task = new Task(QSPRConstants.CROSSVALIDATION, cvConf, getTeacherData(10, 3, null, 0, false, false));
			server.calculateWrapper(task);
			server.allowLocalCalculations = false;
			task.check();
		} finally {
			TransportFactory.waitingTimeBetweenRequests = saveInterval;
		}
	}

	@Test(timeout = MULTIPLIER * 10 * 30000)
	@TaskTest(QSPRConstants.ASNN)
	public void asnnTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.ASNN);
		ASNNConfiguration config = new ASNNConfiguration();
		config.ensemble = 10;
		config.optionsNumber = new ArrayList<Integer>(Arrays.asList(new Integer[] { 2, 2, 1 }));
		Task task = new Task(QSPRConstants.ASNN, config, getTeacherData(25, 5, 15, 3, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task(QSPRConstants.ASNN, wnd.ports.get(1).getValue(0, 0), getTeacherData(10, 5, 10, 3, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 150000)
	@TaskTest(DescriptorsConfiguration.SilicosItScaffold)
	public void inSilicoScaffoldTest() throws Exception {
		DataTable dt = dt_ALOGPS.getSlice(0, 1);
		DataTable dtResult = runTest(new DescriptorsSilicosItScaffoldConfiguration(), dt);
		assertEquals(0, dtResult.errorCount());
	}

	@Test(timeout = MULTIPLIER * 150000)
	@TaskTest(DescriptorsConfiguration.JPLOGP)
	public void jplogpTest() throws Exception {
		DataTable dt = dt_ALOGPS.getSlice(0, 1);
		DataTable dtResult = runTest(new DescriptorsJPlogPConfiguration(), dt);
		assertEquals(0, dtResult.errorCount());
	}

	@Test(timeout = MULTIPLIER * 300000)
	@TaskTest(DescriptorsConfiguration.CDDD)
	public void cdddTest() throws Exception {
		DataTable dt = dt_ALOGPS;
		DataTable dtResult = runTest(new DescriptorsCDDDConfiguration(), dt);
		assertEquals(2, dtResult.errorCount());
	}

	@Test(timeout = MULTIPLIER * 300000)
	@TaskTest(DescriptorsConfiguration.PADEL2)
	public void padel2Test() throws Exception {
		DataTable dt = dt_ALOGPS;
		DescriptorsPaDEL2Configuration cfg = new DescriptorsPaDEL2Configuration();
		DataTable dtResult = runTest(cfg, dt);
		assertEquals(0, dtResult.errorCount());
	}

	@Test(timeout = MULTIPLIER * 300000)
	@TaskTest(DescriptorsConfiguration.MORDRED)
	public void mordredTest() throws Exception {
		DataTable dt = dt_ALOGPS;
		DescriptorsMORDREDConfiguration c = new DescriptorsMORDREDConfiguration();
		c.mordred = Descr.D3;
		DataTable dtResult = runTest(c, dt);
		assertEquals(0, dtResult.errorCount());
	}

	@Test(timeout = MULTIPLIER * 100000)
	@TaskTest(QSPRConstants.COMPARE_SCAFFOLDS)
	public void compareScaffoldsTest() throws Exception {
		WorkflowNodeData twoSets = new WorkflowNodeData();
		twoSets.addPort(getTeacherData(21, 10, 10, 0, false, false).ports.get(0));
		twoSets.addPort(getTeacherData(35, 10, 10, 0, false, false).ports.get(0));

		DataTable dtResult = runTest(QSPRConstants.COMPARE_SCAFFOLDS, new CompareScaffoldsConfiguration(), twoSets);
		assertEquals(0, dtResult.errorCount());
	}

	@Test(timeout = MULTIPLIER * 300000)
	@TaskTest(QSPRConstants.DEEP)
	public void deepTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.DEEP);
		DEEPConfiguration config = new DEEPConfiguration();
		config.optionsNumber = new ArrayList<Integer>(Arrays.asList(new Integer[] { 2, 2, 1 }));
		Task task = new Task(QSPRConstants.DEEP, config, getTeacherData(200, 3, 150, 3, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);

		task = new Task(QSPRConstants.DEEP, wnd.ports.get(1).getValue(0, 0), getTeacherData(10, 3, 1, 3, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}


	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(QSPRConstants.KPLS)
	public void kplsTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.KPLS);
		KPLSConfiguration config = new KPLSConfiguration();
		config.latentVariables = 1;
		Task task = new Task(QSPRConstants.KPLS, config, getTeacherData(200, 3, 100, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task(QSPRConstants.KPLS, wnd.ports.get(1).getValue(0, 0), getTeacherData(10, 3, 1, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(QSPRConstants.KNN)
	public void knnTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.KNN);
		KNNConfiguration config = new KNNConfiguration();
		config.knn = 2;
		Task task = new Task(QSPRConstants.KNN, config, getTeacherData(20, 3, 10, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task(QSPRConstants.KNN, wnd.ports.get(1).getValue(0, 0), getTeacherData(1, 3, 1, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(QSPRConstants.PLS)
	public void plsTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.PLS);
		PLSConfiguration plsc = new PLSConfiguration();
		plsc.numLatentVariables = 5;
		Task task = new Task(QSPRConstants.PLS, plsc, getTeacherData(20, 5, 15, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task(QSPRConstants.PLS, wnd.ports.get(1).getValue(0, 0), getTeacherData(1, 5, 1, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		server = ServerPool.getInstance().getFreeServer(QSPRConstants.PLS);
		task = new Task(QSPRConstants.PLS, new PLSConfiguration(), getTeacherData(20, 3, 15, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		wnd = WorkflowNodeData.fromTask(task);
		task = new Task(QSPRConstants.PLS, wnd.ports.get(1).getValue(0, 0), getTeacherData(1, 3, 1, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(QSPRConstants.MLRA)
	public void mlraTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.MLRA);
		MLRAConfiguration config = new MLRAConfiguration();
		config.optionsNumber = new ArrayList<Integer>(Arrays.asList(new Integer[] { 1 }));

		Task task = new Task(QSPRConstants.MLRA, config, getTeacherData(25, 3, 10, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task(QSPRConstants.MLRA, wnd.ports.get(1).getValue(0, 0), getTeacherData(1, 3, 1, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(QSPRConstants.RFR)
	public void rTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.RFR);
		RFRConfiguration config = new RFRConfiguration();
		config.optionsNumber = new ArrayList<Integer>(Arrays.asList(new Integer[] { 1 }));

		config.numTrees = 10;
		config.numFeatures = 2;
		config.maxDepth = 10;

		Task task = new Task(QSPRConstants.RFR, config, getTeacherData(15, 3, 10, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task(QSPRConstants.RFR, wnd.ports.get(1).getValue(0, 0), getTeacherData(1, 3, 1, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(QSPRConstants.FSMLR)
	public void fsmlrTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.FSMLR);
		Task task = new Task(QSPRConstants.FSMLR, new FSMLRConfiguration(), getTeacherData(20, 3, 15, 0, false, false));
		server.calculateWrapper(task);
		task.check();

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task(QSPRConstants.FSMLR, wnd.ports.get(1).getValue(0, 0), getTeacherData(1, 3, 1, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(QSPRConstants.Workflow)
	public void workflowTest() throws Exception {
		// This test uses workflow with two "+1" incrementors.
		ServerPool.getInstance().servers.add(new AddServer());
		DataTable dtResult = runTest(QSPRConstants.Workflow, new WorkflowConfiguration("TwoIncrements"),
				new DataTable(Integer.valueOf(1)));
		// 1 + 1 + 1 = 3 ?!
		assertEquals(Integer.valueOf(3), dtResult.getValue());
	}

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(DescriptorsConfiguration.RDKIT)
	public void rdkitnewTest() throws Exception {
		DescriptorsRDKITConfiguration rdkit = new DescriptorsRDKITConfiguration();
		DataTable dtResult = runTest(DescriptorsConfiguration.RDKIT, rdkit, wnd_molecules);
		assertTrue(dtResult.containsColumn("MolMR"));
	}

	/**
	 * Checks that number of descriptors for each block == number of descriptors for
	 * all blocks
	 * 
	 * @param configuration
	 * @throws Exception
	 */
	void testDragonConsistency(DescriptorsAbstractDragonConfiguration configuration) throws Exception {
		testDragonConsistency(configuration.getDefaultTypeName(), configuration);
	}

	void testDragonConsistency(String descriptor, DescriptorsAbstractDragonConfiguration configuration)
			throws Exception {

		for (int trial = 0; trial < 2; trial++) {
			int total = 0;
			WorkflowNodeData molecules = wnd_molecules.getDeeperCopy();

			if (trial == 0) { // check whether we get correct number of descriptors for our proxy molecule
				DataTable mol = molecules.ports.get(0);
				for (int i = 0; i < mol.getRowsSize(); i++)
					mol.setValue(i, 0, "C");
			}

			for (int i = 0; i < configuration.getBlocks(); i++) {
				configuration.dragonBlocks = 0;
				configuration.setBit(i);
				DataTable dtResult = runTest(descriptor, configuration, molecules);
				logger.info(descriptor + " block: " + (i + 1) + " columns: " + dtResult.getColumnsSize());
				total += dtResult.getColumnsSize();
			}

			configuration.dragonBlocks = 0;
			for (int i = 0; i < configuration.getBlocks(); i++)
				configuration.setBit(i);
			DataTable dtResult = runTest(descriptor, configuration, wnd_molecules);
			logger.info(descriptor + " all - columns: " + dtResult.getColumnsSize());
			assertTrue(dtResult.getColumnsSize() == total);
			assertTrue(total == configuration.getAllDescriptorsNumber());
		}

	}


	@Test(timeout = MULTIPLIER*150000)
	@TaskTest(DescriptorsConfiguration.CDK2) 
	public void cdk2Test() throws Exception { 
		DataTable dtResult = runTest(new DescriptorsCDK2Configuration(), wnd_molecules); 
		assertTrue(dtResult.containsColumn("ALogP"));
		assertEquals(4.0, dtResult.getValue("Kier1")); 
	}

	@Test(timeout = MULTIPLIER * 150000)
	@TaskTest(DescriptorsConfiguration.MolPrint)
	public void molPrintTest() throws Exception {
		DescriptorsMolPrintConfiguration m = new DescriptorsMolPrintConfiguration();
		m.depth = 3;
		DataTable dtResult = runTest(new DescriptorsMolPrintConfiguration(), "dt_test_MolPrint.xml");
		System.out.println("columns =" + dtResult.getColumnsSize());
		assertTrue(dtResult.getColumnsSize() == 7);
	}

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(QSPRConstants.J48)
	public void j48Test() throws Exception {
		J48Configuration j48conf = new J48Configuration();
		basicWekaTest(QSPRConstants.J48, j48conf, true);
	}

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(QSPRConstants.RF)
	public void rfTest() throws Exception {
		RFConfiguration rfconf = new RFConfiguration();
		rfconf.optionsNumber = new ArrayList<Integer>();
		rfconf.optionsNumber.add(2);

		basicWekaTest(QSPRConstants.RF, rfconf, false);
	}

	@Test(timeout = MULTIPLIER * 300000)
	@TaskTest(QSPRConstants.LIBSVM)
	public void libSvmTest() throws Exception {
		LibSvmConfiguration lsvmconf = new LibSvmConfiguration();
		//
		lsvmconf.gridSearch = true;
		lsvmconf.useWeighting = true;
		lsvmconf.costMax = 1d;
		lsvmconf.costMin = 0d;
		lsvmconf.costStep = 1d;
		lsvmconf.gammaMax = 1d;
		lsvmconf.gammaMin = 1d;
		lsvmconf.gammaStep = 1d;
		lsvmconf.classWeightRatioMax = 0.1;
		lsvmconf.classWeightRatioMin = 0.1;
		lsvmconf.classWeightRatioStep = 0.1;
		lsvmconf.svrEpsilonMax = 1d;
		lsvmconf.svrEpsilonMin = 1d;
		lsvmconf.svrEpsilonStep = 1d;
		//
		lsvmconf.oneClass = false;
		lsvmconf.oneClassLabel = "1";
		lsvmconf.optionsNumber = new ArrayList<Integer>();
		lsvmconf.optionsNumber.add(2);
		lsvmconf.scaleTypeX = ScalingType.RANGE;

		WorkflowNodeData teacherData = getTeacherData(50, 10, 40, 1, true, false);
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.LIBSVM);

		Task task = new Task(QSPRConstants.LIBSVM, lsvmconf, teacherData);
		server.calculateWrapper(task);
		task.check();
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		assertEquals(0, wnd.ports.get(0).errorCount());

		// apply model
		teacherData.ports.remove(1);
		task = new Task(QSPRConstants.LIBSVM, wnd.ports.get(1).getValue(0, 0), teacherData);
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 300000)
	@TaskTest(QSPRConstants.XGBOOST)
	public void gboostTest() throws Exception {
		XGBOOSTConfiguration conf = new XGBOOSTConfiguration();
		conf.depth = 2;
		conf.eta = 1.f;
		conf.rounds = 2;
		conf.objective = "binary:logistic";

		conf.optionsNumber = new ArrayList<Integer>();

		WorkflowNodeData teacherData = getTeacherData(50, 10, 40, 1, true, false);
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.XGBOOST);

		Task task = new Task(QSPRConstants.XGBOOST, conf, teacherData);
		server.calculateWrapper(task);
		task.check();
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		assertEquals(0, wnd.ports.get(0).errorCount());

		// apply model
		teacherData.ports.remove(1);
		task = new Task(QSPRConstants.XGBOOST, wnd.ports.get(1).getValue(0, 0), teacherData);
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 12000000)
	@TaskTest(QSPRConstants.LSSVMG)
	public void lssvmgTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.LSSVMG);
		LSSVMGConfiguration config = new LSSVMGConfiguration();

		config.kernel = "rbf";
		config.cv = 5;

		Task task = new Task(QSPRConstants.LSSVMG, config, getTeacherData(20, 3, 10, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task(QSPRConstants.LSSVMG, wnd.ports.get(1).getValue(0, 0),
				getTeacherData(1, 3, 1, 0, false, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 600000)
	@TaskTest(QSPRConstants.DNN)
	public void dnnTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.DNN);
		DNNConfiguration config = new DNNConfiguration();

		config.epochs = 5;
		config.batchSize = 40;
		// config.modeltype = "dense7new";
		config.modeltype = "dense_exp";

		config.optionsNumber = new ArrayList<Integer>(Arrays.asList(new Integer[] { 2, 2, 1 }));
		Task task = new Task(QSPRConstants.DNN, config, getTeacherData(30, 5, 15, 3, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task(QSPRConstants.DNN, wnd.ports.get(1).getValue(0, 0), getTeacherData(10, 5, 10, 3, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 150000)
	@TaskTest(DescriptorsConfiguration.PyDescriptor)
	public void pyDescriptorTest() throws Exception {
		DataTable dtResult = runTest(new DescriptorsPyDescriptorConfiguration(), dt_PYMOL);
		assertEquals(0, dtResult.errorCount());
		// assertEquals(2, dtResult.getColumnsSize());
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(DescriptorsConfiguration.FRAGMENTS)
	public void fragmentorTest() throws Exception {
		DataTable dtResult = runTest(new DescriptorsFragmentorConfiguration(), "dt_molecules_fragmentor.xml");
		assertTrue(dtResult.containsColumn("C-C"));
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(DescriptorsConfiguration.Mera)
	public void meraTest() throws Exception {
		DataTable dtResult = runTest(new DescriptorsMeraConfiguration(), "dt_molecules_mera.xml");
		assertTrue(dtResult.containsColumn("VME")); // Lame test... make a better one
	}

	@Test(timeout = MULTIPLIER * 1000)
	@TaskTest(DescriptorsConfiguration.RANDOM)
	public void randomTest() {
		// do nothing, just a stub
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(DescriptorsConfiguration.Mersy)
	public void mersyTest() throws Exception {
		DataTable dtResult = runTest(new DescriptorsMersyConfiguration(), "dt_molecules_mera.xml");
		assertTrue(dtResult.containsColumn("SYMC2X")); // Lame test... make a better one
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(DescriptorsConfiguration.GSFrag)
	public void GSFragTest() throws Exception {
		DescriptorsGSFragConfiguration configuration = new DescriptorsGSFragConfiguration();
		configuration.gsfrag = true;
		configuration.gsfragl = true;

		DataTable dtResult = runTest(configuration, wnd_molecules);
		assertTrue(dtResult.containsColumn("p1"));
	}

	@Test(timeout = MULTIPLIER * 150000)
	@TaskTest(DescriptorsConfiguration.ECFP)
	public void ecfpTest() throws Exception {
		DescriptorsECFPConfiguration conf = new DescriptorsECFPConfiguration();
		conf.MORGAN_NBITS = 1024;
		conf.MORGAN_RADIUS = 2;
		DataTable dtResult = runTest(conf, wnd_molecules);
		assertEquals(1024, dtResult.getColumnsSize());
		assertTrue(!dtResult.getRow(0).isError());
	}

	@Test(timeout = MULTIPLIER * 120000)
	@TaskTest(DescriptorsConfiguration.MOPAC)
	public void mopacTestBasicDescriptors() throws Exception {
		DescriptorsMOPACConfiguration conf = new DescriptorsMOPACConfiguration();
		DataTable res = runTest(conf, dt_mopac);
		assertEquals(25, res.getColumnsSize());
		assertEquals(true, res.getRowsNoErrorsSize() == 1);
	}

	@Test(timeout = MULTIPLIER * 120000)
	@TaskTest(DescriptorsConfiguration.MOPAC2016)
	public void mopac2016TestBasicDescriptors() throws Exception {
		DescriptorsMOPAC2016Configuration conf = new DescriptorsMOPAC2016Configuration();
		DataTable res = runTest(conf, dt_mopac);
		assertEquals(35, res.getColumnsSize());
		assertEquals(true, res.getRowsNoErrorsSize() == 1);
	}

	@Test(timeout = MULTIPLIER * 120000)
	@TaskTest(DescriptorsConfiguration.SIGMA)
	public void sigmaBasicDescriptors() throws Exception {
		DescriptorsSIGMAConfiguration conf = new DescriptorsSIGMAConfiguration();
		conf.all = true;
		DataTable res = runTest(conf, dt_mopac);
		assertEquals(163, res.getColumnsSize());
		assertEquals(true, res.getRowsNoErrorsSize() == 1);
	}

	@Test(timeout = MULTIPLIER * 120000)
	@TaskTest(DescriptorsConfiguration.KRAKENX)
	public void krakenxTestBasicDescriptors() throws Exception {
		DescriptorsKrakenXConfiguration conf = new DescriptorsKrakenXConfiguration();
		DataTable res = runTest(conf, dt_mopac);
		assertEquals(124, res.getColumnsSize());
		assertEquals(true, res.getRowsNoErrorsSize() == 1);
	}

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(QSPRConstants.CONSENSUS)
	public void consensusTest() throws Exception {
		assertTrue(true);
	}

	// ///////////////////////////
	// Inductive descriptors //
	// ///////////////////////////

	private void inductiveDescriptorsTest(final String filename, final float[] exp) throws Exception {
		final DataTable res = runTest(DescriptorsConfiguration.INDUCTIVE,
				new DescriptorsInductiveDescriptorsConfiguration(), filename);

		final float[] values = res.asFloatArray()[0];
		final boolean[] pass = new boolean[values.length];
		for (int i = 0; i < values.length; ++i) {
			pass[i] = exp[i] != exp[i] ? true : (Math.abs(values[i] - exp[i]) < 0.001); // NaN means pass because there
			// is no test value available.
			if (!pass[i]) {
				logger.info(
						String.format("'Inductive' descriptor test failure for descriptor %d (expected %f, was %f).",
								i + 1, exp[i], values[i]));
				assertTrue(false);
			}
		}

		assertEquals(54, res.getColumnsSize());
	}

	@Test(timeout = MULTIPLIER * 6000)
	@TaskTest(DescriptorsConfiguration.INDUCTIVE)
	public void inductiveDescriptorsTest() throws Exception {
		inductiveDescriptorsTest("dt_test_inductive_1.xml", new float[] { 2.4735956f, 2.2260001f, 3.2533333f, // Electronegativity
				0.24677353f, 46.904678f, 3.1377628f, 0.89607668f, 0.018982578f, 0.31377628f, 0.29869223f, 0.2063005f,
				0.2756719f, 0.48469195f, 0.31116331f, 0.48469195f, 0.2756719f, // Hardness
				47.119316f, 37.044353f, 10.074966f, 3.624563f, 3.7044351f, 3.3583219f, 2.0631661f, 3.2137465f,
				4.8472981f, 3.6275008f, 2.0631661f, 3.6275008f, // Softness
				1.5245433f, 0, 0.076227166f, -0.25409055f, 0.13188316f, -0.28837961f, // Charge
				2.2148435f, 4.8268681f, 0.64218998f, -0.57556206f, 1.1010636f, -0.19317733f, 3.5208559f, -1.3060123f, // Sigma
				1.2145442f, 0.48443082f, Float.NaN, 0.18009861f, Float.NaN, Float.NaN, Float.NaN, Float.NaN, 1.2145442f,
				0.48443082f, 1.3897711f, 1.1224425f // Steric
		});

		inductiveDescriptorsTest("dt_test_inductive_2.xml", new float[] { 2.2930908f, 2.1233332f, 2.4214287f, // Electronegativity
				0.20669597f, 63.278969f, 3.9921672f, 1.449874f, 0.012918498f, 0.44357413f, 0.20712487f, 0.20492846f,
				0.16547731f, 0.49581984f, 0.27488089f, 0.49581984f, 0.27488089f, // Hardness
				56.255569f, 21.793669f, 34.461903f, 3.5159731f, 2.4215186f, 4.9231291f, 2.0168617f, 3.6379395f,
				4.8797517f, 6.0431252f, 2.0168617f, 3.6379395f, // Softness
				0.81647092f, 0.f, 0.045359496f, -0.058319353f, 0.11969734f, -0.30830145f, // Charge
				1.3500749f, 2.6849072f, 0.57714379f, -0.61248016f, 1.1799701f, -0.11896112f, 2.0174911f, -0.66741616f, // Sigma
				1.0902656f, 0.65316254f, Float.NaN, 0.19258992f, Float.NaN, Float.NaN, Float.NaN, Float.NaN, 1.0401784f,
				0.65316254f, 1.3996975f, 1.1658071f // Steric
		});

		inductiveDescriptorsTest("dt_test_inductive_3.xml", new float[] { 2.3935652f, 2.2145455f, 3.165f, // Electronegativity
				0.25535962f, 49.483391f, 3.6852207f, 0.57038993f, 0.019643048f, 0.33502007f, 0.28519496f, 0.20906724f,
				0.27917686f, 0.49517927f, 0.29121307f, 0.49517927f, 0.27917686f, // Hardness
				45.535f, 38.519127f, 7.0158706f, 3.5026922f, 3.5017388f, 3.5079353f, 2.0194707f, 3.4339118f, 4.7831502f,
				3.581959f, 2.0194707f, 3.581959f, // Softness
				1.0933117f, 0f, 0.049695987f, -0.27332792f, 0.12508838f, -0.29507229f, // Charge
				1.9499263f, 3.6937079f, 0.60668713f, -0.58810419f, 1.1789086f, -0.15117545f, 2.8218172f, -0.87189078f, // Sigma
				1.1324857f, 0.35144436f, Float.NaN, 0.18035793f, Float.NaN, Float.NaN, Float.NaN, Float.NaN, 1.1324857f,
				0.35144436f, 1.3957614f, 1.1280624f // Steric
		});

	}

	@Test(timeout = MULTIPLIER * 12000)
	@TaskTest(QSPRConstants.LEVERAGE)
	public void leverageTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.LEVERAGE);
		DataTable dtDescriptors = getTeacherData(120, 20, 120, 1, false, false).ports.get(0);
		Task task = new Task(QSPRConstants.LEVERAGE, null, new WorkflowNodeData(dtDescriptors), false); // the test
		// modifies
		// dtDescriptors
		// and thus we
		// should create
		// a copy in
		// MongoDB
		server.calculateWrapper(task);
		task.check();
		WorkflowNodeData wndResult = WorkflowNodeData.fromTask(task);

		Serializable conf = wndResult.ports.get(1).getValue(0, 0);

		Double val = (Double) wndResult.ports.get(0).getValue(0, 0);
		task = new Task(QSPRConstants.LEVERAGE, conf, new WorkflowNodeData(dtDescriptors));
		server.calculateWrapper(task);
		task.check();
		wndResult = WorkflowNodeData.fromTask(task);

		assertEquals(val, wndResult.ports.get(0).getValue(0, 0));
	}


	@Test(timeout = MULTIPLIER * 10 * 120000)
	@TaskTest(QSPRConstants.EAGCNG)
	public void eagcngMultiTest() throws Exception {
		String taskname = QSPRConstants.EAGCNG;
		EAGCNGConfiguration conf = new EAGCNGConfiguration();
		conf.nepochs = 10;
		runSMILESAlsoAnalysis(conf, taskname, true, null);
	}

	@Test(timeout = MULTIPLIER * 10 * 120000)
	@TaskTest(QSPRConstants.CHEMPROP)
	public void chempropRegrTest() throws Exception {
		String taskname = QSPRConstants.CHEMPROP;
		ChemPropConfiguration conf = new ChemPropConfiguration();
		conf.nepochs = 2;
		runSMILESAlsoAnalysis(conf, taskname, false, null);
	}

	@Test(timeout = MULTIPLIER * 10 * 120000)
	@TaskTest(QSPRConstants.CHEMPROP)
	public void chempropClassTest() throws Exception {
		String taskname = QSPRConstants.CHEMPROP;
		ChemPropConfiguration conf = new ChemPropConfiguration();
		conf.nepochs = 2;
		runSMILESAnalysis(conf, taskname, true, true, null);
	}

	@Test(timeout = MULTIPLIER * 10 * 120000)
	@TaskTest(QSPRConstants.GNN)
	public void gnnRegrTest() throws Exception {
		String taskname = QSPRConstants.GNN;
		GNNConfiguration conf = new GNNConfiguration();
		conf.nepochs = 2;
		conf.batch = 1024;
		runSMILESAlsoAnalysis(conf, taskname, false, null);
	}

	@Test(timeout = MULTIPLIER * 10 * 120000)
	@TaskTest(QSPRConstants.GNN)
	public void ginClassTest() throws Exception {
		String taskname = QSPRConstants.GNN;
		GNNConfiguration conf = new GNNConfiguration();
		conf.nepochs = 2;
		conf.batch = 1024;
		runSMILESAnalysis(conf, taskname, true, true, null);
	}


	@Test(timeout = MULTIPLIER * 10 * 320000)
	@TaskTest(QSPRConstants.HAMNET)
	public void hamnetClassTest() throws Exception {
		String taskname = QSPRConstants.HAMNET;
		HAMNETConfiguration conf = new HAMNETConfiguration();
		conf.nepochs = 2;
		conf.chirality = true;
		runSMILESAlsoAnalysis(conf, taskname, true, null);
	}


	@Test(timeout = MULTIPLIER * 10 * 320000)
	@TaskTest(QSPRConstants.TRANSNN)
	public void transnnRegrTest() throws Exception {
		String taskname = QSPRConstants.TRANSNN;
		TRANSNNConfiguration conf = new TRANSNNConfiguration();
		conf.nepochs = 2;
		conf.batch = 1024;
		runSMILESAlsoAnalysis(conf, taskname, false, null);
	}
	/*
	 * @Test(timeout = MULTIPLIER*10*320000)
	 * 
	 * @TaskTest(QSPRConstants.TRANSNN) public void transnnmultiRegrTest() throws
	 * Exception { String taskname = QSPRConstants.TRANSNN; TRANSNNConfiguration
	 * conf = new TRANSNNConfiguration(); conf.nepochs = 2; conf.batch = 1024;
	 * conf.chirality = false; runSMILESAlsoAnalysis(conf, taskname, true, null); }
	 * 
	 * 
	 * @Test(timeout = MULTIPLIER*10*320000)
	 * 
	 * @TaskTest(QSPRConstants.TRANSNN) public void transnnClassTest() throws
	 * Exception { String taskname = QSPRConstants.TRANSNN; TRANSNNConfiguration
	 * conf = new TRANSNNConfiguration(); conf.nepochs = 2; conf.batch = 1024;
	 * conf.chirality = false; runSMILESAlsoAnalysis(conf, taskname, true, true,
	 * false); }
	 */

	@Test(timeout = MULTIPLIER * 10 * 120000)
	@TaskTest(QSPRConstants.CNF)
	public void cnfMultiTest() throws Exception {
		String taskname = QSPRConstants.CNF;
		CNFConfiguration conf = new CNFConfiguration();
		// conf.setTaxonomy(new SupportsTaxonomy.Taxonomy("Result0,Result1\n",null));
		conf.implicitValues = new Double[] { 0., 0. };
		conf.nepochs = 5;
		//conf.filterSize = 2;
		conf.setAugmentations(-1, 5,false);
		runSMILESAlsoAnalysis(conf, taskname, true, 5);
	}

	@Test(timeout = MULTIPLIER * 10 * 120000)
	@TaskTest(QSPRConstants.IMGQSAR)
	public void imgMultiTest() throws Exception {
		String taskname = QSPRConstants.IMGQSAR;
		IMGQSARConfiguration conf = new IMGQSARConfiguration();
		conf.implicitValues = new Double[] { 0., 0. };
		conf.nepochs = 2;
		conf.batch_size = 2;
		runSMILESAlsoAnalysis(conf, taskname, true, 5);
	}

	@Test(timeout = MULTIPLIER * 10 * 120000)
	@TaskTest(QSPRConstants.ATTFP)
	public void attfpSingleTest() throws Exception {
		String taskname = QSPRConstants.ATTFP;
		ATTFPConfiguration conf = new ATTFPConfiguration();
		conf.nepochs = 5;
		conf.setAugmentations(5, 5,false);
		conf.lngru = false;
		runSMILESAnalysis(conf, taskname, false, true, 0);
	}

	@Test(timeout = MULTIPLIER * 10 * 120000)
	@TaskTest(QSPRConstants.ATTFP)
	public void attfpMultiTest() throws Exception {
		String taskname = QSPRConstants.ATTFP;
		ATTFPConfiguration conf = new ATTFPConfiguration();
		conf.nepochs = 5;
		conf.setAugmentations(5, 5,false);
		runSMILESAlsoAnalysis(conf, taskname, true, 0);
	}


	@Test(timeout = MULTIPLIER * 100 * 120000)
	@TaskTest(QSPRConstants.KGCNN)
	public void kgcnnTest() throws Exception {

		String taskname = QSPRConstants.KGCNN;

		for(KGCNN method:KGCNN.values()){
			if(method != KGCNN.AttFP)continue; // testing only one method to be faster
			KGCNNConfiguration conf = new KGCNNConfiguration();
			conf.method = method;
			conf.nepochs = 5;
			runSMILESAlsoAnalysis(conf, taskname, false, 0);
		}
	}


	@Test(timeout = MULTIPLIER * 10 * 120000)
	@TaskTest(QSPRConstants.DIMENET)
	public void dimenetSingleTest() throws Exception {
		String taskname = QSPRConstants.DIMENET;
		DIMENETConfiguration conf = new DIMENETConfiguration();
		conf.nbepochs = 5;
		conf.setAugmentations(5, 5,false);
		runSMILESAlsoAnalysis(conf, taskname, false, 0);
	}

	@Test(timeout = MULTIPLIER * 1200000 * 10)
	@TaskTest(QSPRConstants.MACAU)
	public void macauCHEMMulti() throws Exception {
		MACAUConfiguration config = new MACAUConfiguration();
		config.adaptive = true;
		config.burnin = 3;
		config.samples = 10;
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.MACAU);

		config.optionsNumber = new ArrayList<Integer>(Arrays.asList(new Integer[] { 2, 2, 1 }));
		Task task = new Task(QSPRConstants.MACAU, config, getTeacherData(25, 5, 15, 3, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task(QSPRConstants.MACAU, wnd.ports.get(1).getValue(0, 0),
				getTeacherData(10, 5, 10, 3, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}


	@Test(timeout = MULTIPLIER * 1200000 * 10)
	@TaskTest(QSPRConstants.DEEPCHEM)
	public void deepCHEMMulti() throws Exception {
		//for (DeepChemMethod method : DeepChemMethod.values()) 
		{

			DeepChemMethod method = DeepChemMethod.TEXTCNN;
			String taskname = QSPRConstants.DEEPCHEM;
			DeepChemConfiguration conf = new DeepChemConfiguration();
			conf.nepochs = 5;
			conf.method = method;
			if (conf.isTEXTCNN())
				conf.setAugmentations(10, 5, false);

			if (conf.isSupportRegression())
				runSMILESAlsoAnalysis(conf, taskname, true, null);
		}
	}

	@Test(timeout = MULTIPLIER * 1200000)
	@TaskTest(QSPRConstants.DLCA)
	public void dlcaCHEMMulti() throws Exception {
		String taskname = QSPRConstants.DLCA;
		DLCAConfiguration conf = new DLCAConfiguration();
		conf.epochs = 2;
		conf.saveModels = false;
		runSMILESAlsoAnalysis(conf, taskname, true, null);
	}

	@Test(timeout = MULTIPLIER * 15000 * 10 * 10)
	@TaskTest(QSPRConstants.SKLEARN)
	public void chainerSKLEARNSingle() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(QSPRConstants.SKLEARN);

		for (SKLMethod method : SKLMethod.values()) {

			if (method == SKLMethod.EXTRA_TREES)
				continue;
			if (method == SKLMethod.RANDOM_FOREST)
				continue;
			if (method == SKLMethod.RANSAC_TREES)
				continue;
			if (method == SKLMethod.RIDGE_CV)
				continue;
			if (method == SKLMethod.ALL_CV)
				continue;
			if (method == SKLMethod.ALL)
				continue;

			SKLConfiguration conf = new SKLConfiguration();
			conf.method = method;
			conf.optionsNumber = new ArrayList<Integer>();

			WorkflowNodeData teacherData = getTeacherData(50, 10, 40, 1, true, false);

			Task task = new Task(QSPRConstants.SKLEARN, conf, teacherData);
			server.calculateWrapper(task);
			task.check();
			WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
			assertEquals(0, wnd.ports.get(0).errorCount());

			// apply model
			teacherData.ports.remove(1);
			task = new Task(QSPRConstants.SKLEARN, wnd.ports.get(1).getValue(0, 0), teacherData);
			server.calculateWrapper(task);
			task.check();
			assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
			break;
		}
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(QSPRConstants.MMPFrag)
	public void mmpFragServerTest() throws Exception {

		DataTable dtResult = runTest(new MMPFragConfiguration(), wnd_molecules);
		assertTrue(dtResult.errorCount() == 0);
		assertEquals(1, dtResult.getRowsSize());
	}

	// Common utilities

	private static DataTable generateDescriptorsDatatable(int rows, int columns) {
		DataTable dt = new DataTable(true);

		for (int col = 0; col < columns; col++)
			dt.addColumn("Descriptor" + (columns == 1 ? "" : col));

		for (int row = 0; row < rows; row++) {
			dt.addRow();
			dt.getCurrentRow().addAttachment(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM, row);
			dt.getCurrentRow().addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, row);

			for (int col = 0; col < columns; col++)
				dt.setValue(col, Math.random());
		}

		return dt;
	}

	private static DataTable generateLabelsDatatable(int rows) {
		DataTable dt = new DataTable(true);
		dt.addColumn("VALUE");
		for (int row = 0; row < rows; row++) {
			dt.addRow();
			dt.setValue(0, Math.random());
		}
		return dt;
	}

	void basicWekaTest(String taskType, AbstractWekaConfiguration configuration, boolean useCostMatrix)
			throws Exception {
		configuration.optionsNumber = new ArrayList<Integer>();
		configuration.optionsNumber.add(2);

		if (useCostMatrix) {
			// Cost matrix testing
			PropertyWeighting pw = new PropertyWeighting("test_property", 1.0);
			pw.classesWeights = new ArrayList<LabelWeighting.ClassWeighting>();
			pw.classesWeights.add(new ClassWeighting("low", 1.0, new Double[] { 0.0, 1.0 }));
			pw.classesWeights.add(new ClassWeighting("high", 1.0, new Double[] { 1.0, 0.0 }));
			configuration.labelWeighting = new LabelWeighting();
			configuration.labelWeighting.addNewProperty(pw);
			configuration.labelWeighting.useCostMatrix = true;
		}

		CalculationServer server = ServerPool.getInstance().getFreeServer(taskType);

		Task task = new Task(taskType, configuration, getTeacherData(20, 2, 10, 2, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task(taskType, wnd.ports.get(1).getValue(0, 0), getTeacherData(1, 2, 1, 2, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	void runSMILESAlsoAnalysis(ModelAbstractConfiguration conf, String taskname, boolean multi, Integer descriptors)
			throws Exception {
		runSMILESAnalysis(conf, taskname, multi, false, descriptors);
	}

	DataTable getData() throws Exception {
		String vals[] = new String[] { // in total 11 +11 + 2 new molecules
				/*0*/			/* formic acid */" 1\n  -OEChem-01111005042D\n\n  5  4  0     0  0  0  0  0  0999 V2000\n    3.7321    0.5600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.5600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660   -0.5600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  3  1  0  0  0  0\n  1  5  1  0  0  0  0\n  2  3  2  0  0  0  0\n  3  4  1  0  0  0  0\nM  END",
				/* error */"1\n2\n3\n 14 13  0  0  0  0            999 V2000\n    0.3282   -0.3560    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0\n   -0.3863    0.0565    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n    1.0427   -0.7685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7572   -0.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.4716   -0.7685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1861   -0.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.9006   -0.7685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6150   -0.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.3295   -0.7685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0440   -0.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.7585   -0.7685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.4729   -0.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0843   -1.0704    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7407    0.3585    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  3  4  1  0  0  0  0\n  4  5  1  0  0  0  0\n  5  6  1  0  0  0  0\n  6  7  1  0  0  0  0\n  7  8  1  0  0  0  0\n  8  9  1  0  0  0  0\n  9 10  1  0  0  0  0\n 10 11  1  0  0  0  0\n 11 12  1  0  0  0  0\n  1 13  1  0  0  0  0\n  1 14  1  0  0  0  0\nM  CHG  2   1   1   2  -1\nM  END",
				/* acetic acid */" 2\n  -OEChem-01111005082D\n\n  8  7  0     0  0  0  0  0  0999 V2000\n    3.7321    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3100    1.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4631    1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6900    0.2131    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  4  1  0  0  0  0\n  1  8  1  0  0  0  0\n  2  4  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  6  1  0  0  0  0\n  3  7  1  0  0  0  0\nM  END",
				/* acrylic acid */" 3\n  -OEChem-01111005522D\n\n  9  8  0     0  0  0  0  0  0999 V2000\n    4.5981    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    1.3700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4631    0.5600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000   -0.3700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  4  1  0  0  0  0\n  1  9  1  0  0  0  0\n  2  4  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  2  0  0  0  0\n  3  6  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5  8  1  0  0  0  0\nM  END",
				/* propionic acid */" 4\n  -OEChem-01111005122D\n\n 11 10  0     0  0  0  0  0  0999 V2000\n    4.5981    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2646    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.4675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6900    0.7869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4631   -0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3100   -0.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  5  1  0  0  0  0\n  1 11  1  0  0  0  0\n  2  5  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  6  1  0  0  0  0\n  3  7  1  0  0  0  0\n  4  8  1  0  0  0  0\n  4  9  1  0  0  0  0\n  4 10  1  0  0  0  0\nM  END",
				/* crotonic acid */" 6\n  -OEChem-01111005592D\n\n 12 11  0     0  0  0  0  0  0999 V2000\n    5.4641   -0.4400    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5981    1.0600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000   -0.4400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321   -0.4400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5981    0.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.6800    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6900    0.0969    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4631   -0.7500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3100   -0.9769    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321   -1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010   -0.1300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  6  1  0  0  0  0\n  1 12  1  0  0  0  0\n  2  6  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  2  0  0  0  0\n  3  7  1  0  0  0  0\n  4  8  1  0  0  0  0\n  4  9  1  0  0  0  0\n  4 10  1  0  0  0  0\n  5  6  1  0  0  0  0\n  5 11  1  0  0  0  0\nM  END",
				/* butyric acid */" 5\n  -OEChem-01111005542D\n\n 14 13  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3110    0.2131    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.5380    1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6910    1.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  6  1  0  0  0  0\n  1 14  1  0  0  0  0\n  2  6  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  7  1  0  0  0  0\n  3  8  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4  9  1  0  0  0  0\n  4 10  1  0  0  0  0\n  5 11  1  0  0  0  0\n  5 12  1  0  0  0  0\n  5 13  1  0  0  0  0\nM  END",
				/* valeric acid */" 7\n  -OEChem-01111006032D\n\n 17 16  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.8671    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3996    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6025    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.5571   -0.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.4040   -0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.1771    0.7869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  7  1  0  0  0  0\n  1 17  1  0  0  0  0\n  2  7  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  8  1  0  0  0  0\n  3  9  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4 10  1  0  0  0  0\n  4 11  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5 12  1  0  0  0  0\n  5 13  1  0  0  0  0\n  6 14  1  0  0  0  0\n  6 15  1  0  0  0  0\n  6 16  1  0  0  0  0\nM  END",
				/* caproic acid */" 8\n  -OEChem-01111006102D\n\n 20 19  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.8671    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.7331    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3996    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6025    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.4685   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.2656   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.0431    0.2131    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.2700    1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.4231    1.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  8  1  0  0  0  0\n  1 20  1  0  0  0  0\n  2  8  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  9  1  0  0  0  0\n  3 10  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4 11  1  0  0  0  0\n  4 12  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5 13  1  0  0  0  0\n  5 14  1  0  0  0  0\n  6  8  1  0  0  0  0\n  6 15  1  0  0  0  0\n  6 16  1  0  0  0  0\n  7 17  1  0  0  0  0\n  7 18  1  0  0  0  0\n  7 19  1  0  0  0  0\nM  END",
				/* heptanoic acid */" 9\n  -OEChem-01111006132D\n\n 23 22  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.8671    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.7331    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    8.5991    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3996    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6025    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.4685   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.2656   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.1316    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.3346    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.2891   -0.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    9.1360   -0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.9091    0.7869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  9  1  0  0  0  0\n  1 23  1  0  0  0  0\n  2  9  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3 10  1  0  0  0  0\n  3 11  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4 12  1  0  0  0  0\n  4 13  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5 14  1  0  0  0  0\n  5 15  1  0  0  0  0\n  6  8  1  0  0  0  0\n  6 16  1  0  0  0  0\n  6 17  1  0  0  0  0\n  7  9  1  0  0  0  0\n  7 18  1  0  0  0  0\n  7 19  1  0  0  0  0\n  8 20  1  0  0  0  0\n  8 21  1  0  0  0  0\n  8 22  1  0  0  0  0\nM  END",
				/* octanoic acid */"10\n  -OEChem-01111006172D\n\n 26 25  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.8671    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.7331    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    8.5991    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    9.4651    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.4685   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.2656   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3996    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6025    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.1316    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.3346    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.2006   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.9976   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    9.7751    0.2131    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   10.0021    1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    9.1551    1.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1 10  1  0  0  0  0\n  1 26  1  0  0  0  0\n  2 10  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3 11  1  0  0  0  0\n  3 12  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4 13  1  0  0  0  0\n  4 14  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5 15  1  0  0  0  0\n  5 16  1  0  0  0  0\n  6  8  1  0  0  0  0\n  6 17  1  0  0  0  0\n  6 18  1  0  0  0  0\n  7  9  1  0  0  0  0\n  7 19  1  0  0  0  0\n  7 20  1  0  0  0  0\n  8 10  1  0  0  0  0\n  8 21  1  0  0  0  0\n  8 22  1  0  0  0  0\n  9 23  1  0  0  0  0\n  9 24  1  0  0  0  0\n  9 25  1  0  0  0  0\nM  END",
				/* pelargonic acid */"11\n  -OEChem-01111006342D\n\n 29 28  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.8671    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.7331    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    8.5991    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    9.4651    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   10.3312    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.4685   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.2656   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.1316    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.3346    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3996    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6025    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.2006   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.9976   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    9.8637    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    9.0666    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   10.0212   -0.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   10.8681   -0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   10.6412    0.7869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1 11  1  0  0  0  0\n  1 29  1  0  0  0  0\n  2 11  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3 12  1  0  0  0  0\n  3 13  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4 14  1  0  0  0  0\n  4 15  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5 16  1  0  0  0  0\n  5 17  1  0  0  0  0\n  6  8  1  0  0  0  0\n  6 18  1  0  0  0  0\n  6 19  1  0  0  0  0\n  7  9  1  0  0  0  0\n  7 20  1  0  0  0  0\n  7 21  1  0  0  0  0\n  8 10  1  0  0  0  0\n  8 22  1  0  0  0  0\n  8 23  1  0  0  0  0\n  9 11  1  0  0  0  0\n  9 24  1  0  0  0  0\n  9 25  1  0  0  0  0\n 10 26  1  0  0  0  0\n 10 27  1  0  0  0  0\n 10 28  1  0  0  0  0\nM  END",
				/*12*/			/* formic acid */" 1\n  -OEChem-01111005042D\n\n  5  4  0     0  0  0  0  0  0999 V2000\n    3.7321    0.5600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.5600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660   -0.5600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.2500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  3  1  0  0  0  0\n  1  5  1  0  0  0  0\n  2  3  2  0  0  0  0\n  3  4  1  0  0  0  0\nM  END",
				/* error */"1\n2\n3\n 14 13  0  0  0  0            999 V2000\n    0.3282   -0.3560    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0\n   -0.3863    0.0565    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n    1.0427   -0.7685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.7572   -0.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.4716   -0.7685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.1861   -0.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.9006   -0.7685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6150   -0.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.3295   -0.7685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0440   -0.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.7585   -0.7685    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.4729   -0.3560    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   -0.0843   -1.0704    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.7407    0.3585    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  3  4  1  0  0  0  0\n  4  5  1  0  0  0  0\n  5  6  1  0  0  0  0\n  6  7  1  0  0  0  0\n  7  8  1  0  0  0  0\n  8  9  1  0  0  0  0\n  9 10  1  0  0  0  0\n 10 11  1  0  0  0  0\n 11 12  1  0  0  0  0\n  1 13  1  0  0  0  0\n  1 14  1  0  0  0  0\nM  CHG  2   1   1   2  -1\nM  END",
				/* acetic acid */" 2\n  -OEChem-01111005082D\n\n  8  7  0     0  0  0  0  0  0999 V2000\n    3.7321    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3100    1.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4631    1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6900    0.2131    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  4  1  0  0  0  0\n  1  8  1  0  0  0  0\n  2  4  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  6  1  0  0  0  0\n  3  7  1  0  0  0  0\nM  END",
				/* acrylic acid */" 3\n  -OEChem-01111005522D\n\n  9  8  0     0  0  0  0  0  0999 V2000\n    4.5981    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    1.3700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4631    0.5600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000   -0.3700    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  4  1  0  0  0  0\n  1  9  1  0  0  0  0\n  2  4  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  2  0  0  0  0\n  3  6  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5  8  1  0  0  0  0\nM  END",
				/* propionic acid */" 4\n  -OEChem-01111005122D\n\n 11 10  0     0  0  0  0  0  0999 V2000\n    4.5981    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.2646    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.4675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6900    0.7869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4631   -0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3100   -0.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  5  1  0  0  0  0\n  1 11  1  0  0  0  0\n  2  5  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  6  1  0  0  0  0\n  3  7  1  0  0  0  0\n  4  8  1  0  0  0  0\n  4  9  1  0  0  0  0\n  4 10  1  0  0  0  0\nM  END",
				/* crotonic acid */" 6\n  -OEChem-01111005592D\n\n 12 11  0     0  0  0  0  0  0999 V2000\n    5.4641   -0.4400    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5981    1.0600    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000   -0.4400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321   -0.4400    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.5981    0.0600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.8660    0.6800    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.6900    0.0969    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    1.4631   -0.7500    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.3100   -0.9769    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.7321   -1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010   -0.1300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  6  1  0  0  0  0\n  1 12  1  0  0  0  0\n  2  6  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  2  0  0  0  0\n  3  7  1  0  0  0  0\n  4  8  1  0  0  0  0\n  4  9  1  0  0  0  0\n  4 10  1  0  0  0  0\n  5  6  1  0  0  0  0\n  5 11  1  0  0  0  0\nM  END",
				/* butyric acid */" 5\n  -OEChem-01111005542D\n\n 14 13  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3110    0.2131    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.5380    1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6910    1.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  6  1  0  0  0  0\n  1 14  1  0  0  0  0\n  2  6  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  7  1  0  0  0  0\n  3  8  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4  9  1  0  0  0  0\n  4 10  1  0  0  0  0\n  5 11  1  0  0  0  0\n  5 12  1  0  0  0  0\n  5 13  1  0  0  0  0\nM  END",
				/* valeric acid */" 7\n  -OEChem-01111006032D\n\n 17 16  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.8671    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3996    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6025    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.5571   -0.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.4040   -0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.1771    0.7869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  7  1  0  0  0  0\n  1 17  1  0  0  0  0\n  2  7  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  8  1  0  0  0  0\n  3  9  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4 10  1  0  0  0  0\n  4 11  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5 12  1  0  0  0  0\n  5 13  1  0  0  0  0\n  6 14  1  0  0  0  0\n  6 15  1  0  0  0  0\n  6 16  1  0  0  0  0\nM  END",
				/* caproic acid */" 8\n  -OEChem-01111006102D\n\n 20 19  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.8671    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.7331    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3996    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6025    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.4685   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.2656   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.0431    0.2131    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.2700    1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.4231    1.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  8  1  0  0  0  0\n  1 20  1  0  0  0  0\n  2  8  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  9  1  0  0  0  0\n  3 10  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4 11  1  0  0  0  0\n  4 12  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5 13  1  0  0  0  0\n  5 14  1  0  0  0  0\n  6  8  1  0  0  0  0\n  6 15  1  0  0  0  0\n  6 16  1  0  0  0  0\n  7 17  1  0  0  0  0\n  7 18  1  0  0  0  0\n  7 19  1  0  0  0  0\nM  END",
				/* heptanoic acid */" 9\n  -OEChem-01111006132D\n\n 23 22  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.8671    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.7331    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    8.5991    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3996    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6025    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.4685   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.2656   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.1316    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.3346    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.2891   -0.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    9.1360   -0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.9091    0.7869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  9  1  0  0  0  0\n  1 23  1  0  0  0  0\n  2  9  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3 10  1  0  0  0  0\n  3 11  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4 12  1  0  0  0  0\n  4 13  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5 14  1  0  0  0  0\n  5 15  1  0  0  0  0\n  6  8  1  0  0  0  0\n  6 16  1  0  0  0  0\n  6 17  1  0  0  0  0\n  7  9  1  0  0  0  0\n  7 18  1  0  0  0  0\n  7 19  1  0  0  0  0\n  8 20  1  0  0  0  0\n  8 21  1  0  0  0  0\n  8 22  1  0  0  0  0\nM  END",
				/* octanoic acid */"10\n  -OEChem-01111006172D\n\n 26 25  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.8671    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.7331    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    8.5991    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    9.4651    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.4685   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.2656   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3996    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6025    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.1316    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.3346    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.2006   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.9976   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    9.7751    0.2131    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   10.0021    1.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    9.1551    1.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1 10  1  0  0  0  0\n  1 26  1  0  0  0  0\n  2 10  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3 11  1  0  0  0  0\n  3 12  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4 13  1  0  0  0  0\n  4 14  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5 15  1  0  0  0  0\n  5 16  1  0  0  0  0\n  6  8  1  0  0  0  0\n  6 17  1  0  0  0  0\n  6 18  1  0  0  0  0\n  7  9  1  0  0  0  0\n  7 19  1  0  0  0  0\n  7 20  1  0  0  0  0\n  8 10  1  0  0  0  0\n  8 21  1  0  0  0  0\n  8 22  1  0  0  0  0\n  9 23  1  0  0  0  0\n  9 24  1  0  0  0  0\n  9 25  1  0  0  0  0\nM  END",
				/* pelargonic acid */"11\n  -OEChem-01111006342D\n\n 29 28  0     0  0  0  0  0  0999 V2000\n    2.5369    0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030   -0.7500    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n    6.8671    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    7.7331    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.0010    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    8.5991    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    5.1350    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    9.4651    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    4.2690    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n   10.3312    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.4030    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    6.4685   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.2656   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.1316    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    7.3346    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    6.3996    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.6025    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.2006   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    8.9976   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.7365   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    5.5335   -0.2249    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    9.8637    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    9.0666    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    4.6675    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    3.8705    1.2250    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   10.0212   -0.2869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   10.8681   -0.0600    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n   10.6412    0.7869    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n    2.0000    0.4400    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1 11  1  0  0  0  0\n  1 29  1  0  0  0  0\n  2 11  2  0  0  0  0\n  3  4  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3 12  1  0  0  0  0\n  3 13  1  0  0  0  0\n  4  6  1  0  0  0  0\n  4 14  1  0  0  0  0\n  4 15  1  0  0  0  0\n  5  7  1  0  0  0  0\n  5 16  1  0  0  0  0\n  5 17  1  0  0  0  0\n  6  8  1  0  0  0  0\n  6 18  1  0  0  0  0\n  6 19  1  0  0  0  0\n  7  9  1  0  0  0  0\n  7 20  1  0  0  0  0\n  7 21  1  0  0  0  0\n  8 10  1  0  0  0  0\n  8 22  1  0  0  0  0\n  8 23  1  0  0  0  0\n  9 11  1  0  0  0  0\n  9 24  1  0  0  0  0\n  9 25  1  0  0  0  0\n 10 26  1  0  0  0  0\n 10 27  1  0  0  0  0\n 10 28  1  0  0  0  0\nM  END",
		"12\n OpenBabel01021914222D\n\n  7  7  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  6  1  0  0  0  0\n  1  2  2  0  0  0  0\n  2  3  1  0  0  0  0\n  3  4  2  0  0  0  0\n  4  5  1  0  0  0  0\n  5  6  2  0  0  0  0\n  6  7  1  0  0  0  0\nM  END" };

		int err = 2, n = 0;

		DataTable data = new DataTable(true);
		for (String mol : vals) {
			AbstractDataRow row = data.addRow();
			row.addAttachment(QSPRConstants.SMILES_ATTACHMENT, Various.molecule.convertToFormat(mol, QSPRConstants.SMILESH));
			if (err-- == 0)
				row.addAttachment(QSPRConstants.SMILES_ATTACHMENT, "ERROR");
			row.addAttachment(QSPRConstants.INCHIKEYS, Various.molecule.getInChiKey(mol));
			row.addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, n++);
		}

		return data;
	}

	void runSMILESAnalysis(ModelAbstractConfiguration conf, String taskname, boolean multi, boolean classonly,
			Integer descriptors) throws Exception {

		WorkflowNodeData wndInput;

		DataTable data = getData();

		if (descriptors != null && descriptors > 0)
			for (int i = 0; i < data.getRowsSize(); i++) {
				for (int j = 0; j < descriptors; j++)
					data.setValue(i, QSPRConstants.DESCRIPTOR + j, i * descriptors + j + 0.1);
			}

		int number = data.getRowsSize();

		if (classonly && !multi){
			wndInput = getTeacherData(number, 1, number, 1, true, false);
			conf.optionsNumber = new ArrayList<Integer>(Arrays.asList(new Integer[] { 2 }));
			wndInput.ports.set(0, data);
		}
		else
			if (multi) {
				wndInput = getTeacherData(number, 1, number, 2, true, false);
				conf.optionsNumber = new ArrayList<Integer>(Arrays.asList(new Integer[] { 2, 2 }));
				wndInput.ports.set(0, data);
			}
			else {
				DataTable labels = new DataTable(true);
				labels.addColumn(QSPRConstants.VALUE);
				for (int i = 0; i < number; i++) {
					labels.addRow();
					labels.setValue(Various.molecule.getMass((String)data.getRow(i).getAttachment(QSPRConstants.SMILES_ATTACHMENT)));
				}
				conf.optionsNumber = new ArrayList<Integer>(Arrays.asList(new Integer[] { 1 }));
				wndInput = new WorkflowNodeData().addPort(data).addPort(labels).addPort(data);
			}

		CalculationServer server = ServerPool.getInstance().getFreeServer(taskname);

		Task task = new Task(taskname, conf, wndInput);
		server.calculateWrapper(task);
		task.check();

		int errors = descriptors != null && descriptors > 0 ? 1 : 2; // fix?

		DataTable res = ((WorkflowNodeData) task.getResult()).ports.get(0);

		System.out.println("errors: " + res.errorCount() + " noerrors= " + res.getRowsNoErrorsSize() +  " all= " + res.getRowsSize() );
		assertTrue( (res.getRowsNoErrorsSize() + errors) == res.getRowsSize()); // two errors because of compressing

		wndInput.ports.remove(1);

		if(conf.saveModels()) {
			task = new Task(taskname, conf, wndInput);
			server.calculateWrapper(task);
			task.check();
			WorkflowNodeData wf = ((WorkflowNodeData) task.getResult());
			res = ((WorkflowNodeData) task.getResult()).ports.get(0);
			System.out.println("errors: " + res.errorCount());
			// assertEquals(errors, wf.ports.get(0).errorCount()); // one error is removed
			// before calculations
			assertEquals(number, wf.ports.get(0).getRowsSize());
		}

	}

	// NOT USED

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest("ASNNP")
	public void asnnpTest() throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer("ASNN");
		ASNNConfiguration config = new ASNNConfiguration();
		config.ensemble = 10;
		config.optionsNumber = new ArrayList<Integer>(Arrays.asList(new Integer[] { 2, 2, 1 }));
		config.additionalParam = ASNNServer.PARALLEL + "=" + config.ensemble / 2;
		Task task = new Task("ASNN", config, getTeacherData(25, 5, 15, 3, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());

		// apply model
		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		task = new Task("ASNN", wnd.ports.get(1).getValue(0, 0), getTeacherData(10, 5, 10, 3, true, false));
		server.calculateWrapper(task);
		task.check();
		assertEquals(0, ((WorkflowNodeData) task.getResult()).ports.get(0).errorCount());
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(QSPRConstants.OBABEL)
	public void babelTest() throws Exception {
		DataTable dtResult = runTest(QSPRConstants.OBABEL, new OpenBabelConfiguration(true), "dt_conformations.xml");
		babelTest(dtResult);
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(QSPRConstants.CORINA)
	public void corinaTest() throws Exception {
		DataTable dtResult = runTest(QSPRConstants.CORINA, new CorinaConfiguration(true), wnd_molecules);
		assertTrue(!dtResult.getRow(0).isError());
	}

	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest(QSPRConstants.OBGEN)
	public void obgenTest() throws Exception {
		DataTable dtResult = runTest(QSPRConstants.OBGEN, new ObgenConfiguration(true), "dt_conformations.xml"); // remote
		// version
		assertTrue(dtResult.errorCount() == 0);
		dtResult = runTest(QSPRConstants.OBGEN, new ObgenConfiguration(true), "dt_conformations.xml"); // local version
		assertTrue(dtResult.errorCount() == 0);
	}

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(QSPRConstants.BALLOON)
	public void balloonTest() throws Exception {
		DataTable dtResult = runTest(QSPRConstants.BALLOON, new BalloonConfiguration(true), "dt_conformations.xml");
		assertTrue(dtResult.errorCount() == 0);
		dtResult = runTest(QSPRConstants.BALLOON, new BalloonConfiguration(true), "dt_conformations.xml");
		assertTrue(dtResult.errorCount() == 0);
	}

	private void babelTest(DataTable dtResult) {
		String lines[] = dtResult.getRow(0).getValue(0).toString().split("\n");
		boolean finishTag = false;
		for (int i = lines.length - 1; i > 0; i--)
			if (lines[i].equals("M  END")) {
				finishTag = true;
				break;
			}

		assertTrue(lines[3].startsWith("  9  9"));
		assertTrue(finishTag);
	}

	@Test(timeout = MULTIPLIER * 30000)
	@TaskTest(DescriptorsConfiguration.StructuralAlerts)
	public void structuralAlertsTest() throws Exception {
		DataTable dtMolecules = new DataTable();
		dtMolecules.addColumn(QSPRConstants.SDF_COLUMN);
		dtMolecules.addRow();
		dtMolecules.setValue(0, 0, Various.molecule.convertToFormat("CCCCNCCCCCCC=C", QSPRConstants.SDF));

		DescriptorsStructuralAlertsConfiguration conf = new DescriptorsStructuralAlertsConfiguration();
		conf.compactMode = false;
		String alerts[] = {"NOT C=C AND NC","NOT NCCC","CC","CCCCC","C=C AND NC AND CCC","CCF OR NCC","CCF OR NOT CCF","NOT CC OR NOT CCN"};
		conf.setAlerts(Arrays.asList(alerts));

		dtMolecules.columnAttachments.put(QSPRConstants.SDF_COLUMN, map);

		DataTable dtResult = runTest(DescriptorsConfiguration.StructuralAlerts, conf, dtMolecules);

		boolean[] matches = new boolean[] { false, false, true, true, true, true, true, false };
		for (int i = 0; i < matches.length; i++)
			assertEquals(matches[i], ((Double) dtResult.getValue(0, i)) > 0.0);
	}

	private WorkflowNodeData getTeacherData(int samples, int dimensions, Integer traingsetsize, int classes,
			boolean withClassification, boolean conditions) throws Exception {
		DataTable leftSide = new DataTable(true);
		DataTable rightSide = new DataTable(true);
		Random r = new Random(100);

		if (conditions)
			leftSide.addColumns(dimensions + 1); // conditions are present, one condition -- testing
		else
			leftSide.addColumns(dimensions); // data itself

		if (classes > 0) {
			rightSide.addColumn(QSPRConstants.VALUE); // Column with data
			rightSide.addColumn(QSPRConstants.CLASS); // column with type of data
		} else
			rightSide.addColumn(QSPRConstants.VALUE); // Column with data

		// generation of random data
		for (int i = 0; i < samples; i++) {
			// generation of input data
			AbstractDataRow row=leftSide.addRow();
			row.addAttachment(QSPRConstants.INCHIKEYS, "ATUOYWHBWRKTHZ-UHFFFAOYSA-N");
			row.addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, Integer.valueOf(i));
			row.addAttachment(QSPRConstants.RECORD_ID_ATTACHMENT, Long.valueOf(i));
			row.addAttachment(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM, Integer.valueOf(i));
			row.addAttachment(QSPRConstants.SMILES_ATTACHMENT, "CCC");
			double sum = 0;
			for (int k = 0; k < dimensions; k++) {
				double a = 10 * r.nextDouble();
				sum += k * a;
				if (dimensions == 1)
					a = i;
				leftSide.setValue(k, Double.valueOf(a));
			}
			if (conditions) {
				int randint = r.nextInt(3);
				leftSide.setValue(dimensions, Long.valueOf(randint));
			}
			sum += r.nextDouble() / 10; // Noise

			// generation of target values
			if (traingsetsize == null || i < traingsetsize) {
				rightSide.addRow();

				if (withClassification) {
					// if we have multiple properties with classification
					int randint = r.nextInt(classes);
					if (randint == 2) {
						// set value for regression
						// logger.info("class "+randint+" sum"+sum);
						rightSide.setValue(0, Double.valueOf(sum));
					} else {
						// set values for classification
						int options = r.nextInt(2);
						// logger.info("class "+randint+" option"+options);
						rightSide.setValue(0, Double.valueOf(options));
					}
					rightSide.setValue(1, Long.valueOf(randint));
				} else if (classes > 0)
					// if we have multiple properties only with regression problem
				{
					int randint = r.nextInt(classes);
					// logger.info(randint);
					rightSide.setValue(0, Double.valueOf(sum));
					rightSide.setValue(1, Long.valueOf(randint));
				} else
					// if we have only one output value
					rightSide.setValue(Double.valueOf(sum));
			} else
				rightSide.addStubRow();

		}

		return new WorkflowNodeData(leftSide).addPort(rightSide);
	}

	private DataTable runTest(String serverName, String taskType, Serializable configuration, WorkflowNodeData wfnd)
			throws Exception {
		CalculationServer server = ServerPool.getInstance().getFreeServer(serverName);
		if (server == null) {
			logger.info("server " + taskType + " is null");
			throw new IOException("server " + taskType + " is null");
		}
		Task task = new Task(taskType, configuration, wfnd);
		server.calculateWrapper(task);
		task.check();

		DataTable dtResult = ((WorkflowNodeData) task.getResult()).ports.get(0);
		return dtResult;
	}


	private DataTable runTest(String taskType, Serializable configuration, WorkflowNodeData wfnd) throws Exception {
		return runTest(taskType, taskType, configuration, wfnd);
	}

	private DataTable runTest(String taskType, Serializable configuration, DataTable dt) throws Exception {
		dt.columnAttachments.put(QSPRConstants.SDF_COLUMN, map);
		return runTest(taskType, taskType, configuration, new WorkflowNodeData(dt));
	}

	private DataTable runTest(String taskType, Serializable configuration, String fileName) throws Exception {
		DataTable dt = DataTable.fromXml(jc, PATH + fileName);
		dt.columnAttachments.put(QSPRConstants.SDF_COLUMN, map);
		return runTest(taskType, taskType, configuration, new WorkflowNodeData(dt));
	}

	private DataTable runTest(DescriptorsAbstractConfiguration configuration, DataTable dt) throws Exception {
		dt.columnAttachments.put(QSPRConstants.SDF_COLUMN, map);
		return runTest(configuration.getDefaultTypeName(), configuration.getDefaultTypeName(), configuration,
				new WorkflowNodeData(dt));
	}

	private DataTable runTest(DescriptorsAbstractConfiguration configuration, WorkflowNodeData wfnd) throws Exception {
		return runTest(configuration.getDefaultTypeName(), configuration.getDefaultTypeName(), configuration, wfnd);
	}

	private DataTable runTest(DescriptorsAbstractConfiguration configuration, String string) throws Exception {
		return runTest(configuration.getDefaultTypeName(), configuration, string);
	}


	@Test(timeout = MULTIPLIER * 15000)
	@TaskTest("ExpValues")
	public void expValuesTest() throws Exception {
		DataTable dt = new DataTable(Various.molecule.convertToFormat("CC=O", QSPRConstants.SDF));
		dt.getRow(0).addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, 4448);

		DescriptorsExpValuesConfiguration conf = new DescriptorsExpValuesConfiguration();
		conf.multipleExpValues = "USE_MAXIMUM";
		conf.basketId = 7428L;
		conf.properties.add("logPow");

		DataTable dtResult = runTest("ExpValues", conf, new WorkflowNodeData(dt));
		assertTrue(Math.abs((Double) dtResult.getValue(0, 0) + 0.34D) < 1E-3);
	}

	@Test(timeout = MULTIPLIER * 180000)
	@TaskTest("RA")
	public void raTest() throws Exception {

		CalculationServer server = ServerPool.getInstance().getFreeServer("RA");

		String str1 = "-2\n-1.7\n-1\n-0.5\n0\n0.5\n1\n1.7\n2\n";
		String str2 = "-1\n-0.7\n0\n0.5\n1\n1.5\n2\n2.7\n3";
		DataTable dtIn = new DataTable("");
		dtIn.addColumn("s1");
		dtIn.setValue("s1", str1);
		dtIn.addColumn("s2");
		dtIn.setValue("s2", str2);

		WorkflowNodeData wfnd = new WorkflowNodeData(dtIn);

		Task task = new Task("RA", null, wfnd);
		server.calculateWrapper(task);
		task.check();
		DataTable dtResult = ((WorkflowNodeData) task.getResult()).ports.get(0);
		int result = Integer.parseInt(dtResult.getValue().toString().substring(0, 2));
		assertTrue(result < 80 && result > 76);

	}

}

// "+1" server, used in Workflow test
class AddServer extends WorkflowNodeServer {
	public AddServer() {
		supportedTaskType = "Increment";
	}

	@Override
	public WorkflowNodeData calculate(WorkflowNodeData task, Serializable configuration) throws Exception {
		Integer result = Integer.valueOf((Integer) task.ports.get(0).getValue(0, 0) + 1);
		return new WorkflowNodeData(new DataTable(result));
	}
}

class TestModelServ extends WorkflowNodeServer {

	public TestModelServ() {
		supportedTaskType = "TestModel";
	}

	@Override
	public WorkflowNodeData calculate(WorkflowNodeData task, Serializable configuration) throws Exception {
		DataTable dtPredictions = new DataTable(true);
		DataTable dtModel = new DataTable(new ASNNConfiguration());
		dtPredictions.addColumns(1);
		DataTable leftSide = task.ports.get(0);
		leftSide.reset();
		while (leftSide.nextRow()) {
			dtPredictions.addRow();
			dtPredictions.setValue(new Double(0));
		}

		return new WorkflowNodeData(dtPredictions).addPort(dtModel);
	}

}