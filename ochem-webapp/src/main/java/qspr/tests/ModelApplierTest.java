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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.MockTransport;
import qspr.metaserver.transport.TransportFactory;
import qspr.workflow.utils.QSPRConstants
;
import qspr.modelling.applier.ModelApplier;
import qspr.util.MoleculePeer;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

@PeriodicTest(schedule = ScheduleType.DAILY, type="general")
public class ModelApplierTest 
{
	private static transient final Logger logger = LogManager.getLogger(ModelApplierTest.class);

	@Rule
	public TestRule loginRule = new LoginRule("MAT", true);

	@Test
	@NotifyOnFailure(minConsequentFailures = 3)
	public void ModelApplier() throws Exception
	{
		logger.info("Starting ModelApplierTestBasicTest");
		Basket ts = OCHEMTestHelpers.generateRandomBasketNumeric(OCHEMTestHelpers.getRand(), 10);
		Model m = OCHEMTestHelpers.trainAModel("MAT", ts, null, "classpath:tests/ann-cv.ochem.xml", true, true);

		MockTransport predictionMock = new PredictionMockTransport();
		TransportFactory.setThreadTransport(predictionMock);

		ModelApplier applier = new ModelApplier();
		applier.compoundsProvider.basket = OCHEMTestHelpers.generateRandomBasketNumeric(OCHEMTestHelpers.getRand(), 10);
		applyModel(applier, m);

		Thread.sleep(10000);
		
		// 10 new predictions
		Assert.assertEquals(1.0, applier.modelTasks.get(0).wndResult.ports.get(0).getValue(0, 0));
		Assert.assertEquals(0, applier.modelTasks.get(0).wndResult.ports.get(0).errorCount());

		applyModel(applier, m);
		
		Thread.sleep(10000); // to wait for synchronization of cache

		// All predictions should have been cached before
		Assert.assertEquals(10, applier.modelTasks.get(0).getCachedCount());
		Assert.assertEquals(-1L, applier.modelTasks.get(0).pTask.taskId.longValue());
		Assert.assertEquals(0, applier.modelTasks.get(0).wndResult.ports.get(0).errorCount());

		// Add one mol and calculate the cache
		Basket newBasket = new Basket();
		newBasket.addEntry(applier.compoundsProvider.basket.entries.get(0).ep);
		newBasket.addMolecule(MoleculePeer.getMolecule("CCCCFCC"));
		applier = new ModelApplier();
		applier.compoundsProvider.setBasket(newBasket);
		m = Repository.model.getById(m.id);
		applyModel(applier, m);
		Assert.assertEquals(1, applier.modelTasks.get(0).getCachedCount());
		Assert.assertEquals(2, applier.modelTasks.get(0).wndResult.ports.get(0).getRowsSize());

		// Simulate fetching from pending task
		applier = new ModelApplier(applier.modelTasks.get(0).pTask);
		Assert.assertEquals(2, applier.modelTasks.get(0).wndResult.ports.get(0).getRowsSize());

		// Now fail and check if cached predictions are reported
		predictionMock.fail = true;
		applier.compoundsProvider.basket.addMolecule(MoleculePeer.getMolecule("CCF=CCC"));
		applyModel(applier, m);
		Assert.assertEquals(2, applier.modelTasks.get(0).getCachedCount());
		Assert.assertEquals(3, applier.modelTasks.get(0).wndResult.ports.get(0).getRowsSize());

		// Simulate fetching from pending task
		applier = new ModelApplier(applier.modelTasks.get(0).pTask);
		Assert.assertEquals(3, applier.modelTasks.get(0).wndResult.ports.get(0).getRowsSize());

		ModelApplier.clearCachedPredictions(m.publicId);
	}

	/**
	 * Simulate a failure and check that applier reacts appropriately
	 */
	@Test
	public void failuresTest() throws Exception {
		Basket ts = OCHEMTestHelpers.generateRandomBasketNumeric(OCHEMTestHelpers.getRand(), 10);
		Model m = OCHEMTestHelpers.trainAModel("MAT", ts, null, "classpath:tests/ann-cv.ochem.xml", true, true);

		MockTransport predictionMock = new PredictionMockTransport();
		TransportFactory.setThreadTransport(predictionMock);
		predictionMock.fail = true;

		ModelApplier applier = new ModelApplier();
		applier.compoundsProvider.basket.addMolecule("CCC");
		applyModel(applier, m);
		Assert.assertTrue(applier.isError());
	}

	/**
	 * A helper model to configure applier and run the predictions.
	 */
	private void applyModel(ModelApplier applier, Model m) throws Exception
	{
		applier.modelTasks.clear();
		applier.addModel(m);
		applier.modelTasks.get(0).setUseCache();
		Globals.restartAllTransactions(true);
		applier.startAndWait();
	}
}

class PredictionMockTransport extends MockTransport
{
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public Serializable getTaskResult(Task inTask) throws IOException, ClassNotFoundException
	{
		WorkflowNodeData in = (WorkflowNodeData) inTask.getData();
		DataTable dtResult = new DataTable(true);
		dtResult.addColumn(QSPRConstants.PREDICTION_RESULT_COLUMN);
		for (int i = 0; i < in.ports.get(0).getRowsSize(); i++)
			dtResult.addRow().setValue(0, 1.0d);

		return new WorkflowNodeData(dtResult);
	}
}
