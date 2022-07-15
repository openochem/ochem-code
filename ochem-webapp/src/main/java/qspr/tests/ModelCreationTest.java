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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;

import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import qspr.Environment;
import qspr.Globals;
import qspr.ThreadScope;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.Property;
import qspr.export.ExportableSetConfiguration;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.modelling.ModelStatistics;
import qspr.modelling.SetStatistics;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.applier.ModelApplierTaskProcessor;
import qspr.util.ExportThread;
import qspr.util.WrapperThread;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.utils.NumericalValueStandardizer;

/**
 * Test the full modeling workflow.
 * 
 * The test is based on the externally saved test case files:
 * - training set (training-set.csv)
 * - validation set (validation-set.csv)
 * - model configuration file (model.ochem.xml)
 * - expected statistical parameters (statistics.ochem.xml)
 * - expected prediction values for the validation set (predictions.txt)
 * 
 * Expected location of the files: WEB-INF/tests/models and sub-folders for models
 * 
 * @author midnighter
 *
 */
@PeriodicTest(schedule = ScheduleType.TWELVE_HOURS, type = "general")
public class ModelCreationTest 
{
	public static boolean generateReferenceFiles;
	private static transient Logger logger = LogManager.getLogger(ModelCreationTest.class);

	private static final String PREFIX = "MCT";

	/**
	 * Allowed percentage of inconsistent predictions (30%)
	 */
	private static double INCONSISTENCY_TOLERACE = 0.3;

	@BeforeClass
	public static void before()
	{
		logger.info("Deleting the garbage from previous tests...");
		new WrapperThread() {

			@Override
			public void wrapped() throws Exception {
				// Delete garbage from previous tests
				OCHEMTestHelpers.cleanupByPrefix(PREFIX);
			}
		}.run();
	}

	@AfterClass
	public static void after()
	{
		logger.info("Deleting previous tests...");
		new WrapperThread() {

			@Override
			public void wrapped() throws Exception {
				// Delete garbage from previous tests
				OCHEMTestHelpers.cleanupByPrefix(PREFIX);
			}
		}.run();
	}

	@Test(timeout = 600000)
	@NotifyOnFailure(minConsequentFailures = 1)
	public void createAndApplyCV_ASNN_ModelTest() throws Exception
	{
		createAndApplyModelTest("ann-cv-oestate");
	}

	@Test(timeout = 600000)
	@NotifyOnFailure(minConsequentFailures = 3)
	public void createAndApplyCV_SVM_ModelTest() throws Exception
	{
		createAndApplyModelTest("libsvm-cv-oestate");
	}

	@Test(timeout = 900000)
	@NotifyOnFailure(minConsequentFailures = 1)
	public void createAndApplyBagging_ASNN_ModelTest() throws Exception
	{
		createAndApplyModelTest("ann-bagging-oestate");
	}

	public static void createAndApplyModelTest(final String testModelFolder) throws Exception
	{
		createAndApplyModelTest(testModelFolder, Property.TYPE_NUMERIC);
	}

	public static void createAndApplyModelTest(final String testModelFolder, final int propertyType) throws Exception
	{
		ThreadScope.get().threadClassLoader = ModelCreationTest.class.getClassLoader();
		WrapperThread t = new DataDrivenTestWrapper(PREFIX) 
		{
			@Override
			public void wrappedTest() throws Exception 
			{
				logger.info("Starting modeling process for " + testModelFolder);

				String testBasket = getTestsFolder();
				String testModel = getTestsFolder() + "/" + testModelFolder + "/model.ochem.xml";

				Model model = OCHEMTestHelpers.trainAModel(PREFIX, testBasket, testModel, propertyType, false, false);
				ModelStatistics mStats = (ModelStatistics) model.modelMappings.get(0).statisticsOriginal.getObject();

				String statisticsPath = getTestsFolder() + "/" + testModelFolder + "/statistics.ochem.xml";
				if (generateReferenceFiles)
				{
					logger.info("Saving statistics to " + statisticsPath);
					for (SetStatistics ss : mStats.sets)
						ss.points.clear();
					marshal(mStats, statisticsPath);
				}
				else
					compareStatistics(statisticsPath, mStats);

				testApplier(model, testModelFolder);
			}
		};
		t.run();

		if (t.exception != null)
			throw t.exception;
	}

	private static void marshal(Object obj, String fileName) throws JAXBException
	{
		Marshaller marshaller = Globals.jaxbContext.createMarshaller();
		marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, new Boolean(true));
		marshaller.marshal(obj, new File(fileName));
	}

	/**
	 * Compare the statistics saved to a file with the newly calculated statistics
	 */
	private static void compareStatistics(String exampleFile, ModelStatistics msActual) throws JAXBException, IOException
	{
		ModelStatistics msExpected = (ModelStatistics) Globals.jaxbContext.createUnmarshaller().unmarshal(new File(exampleFile));

		for (int i = 0; i < msActual.sets.size(); i++)
		{
			SetStatistics ssExpected = msExpected.sets.get(i);
			SetStatistics ssActual = msActual.sets.get(i);
			Assert.assertEquals(ssExpected.r2, Double.valueOf(ssActual.getR2Str()), ssActual.r2Std); 
			Assert.assertEquals(ssExpected.q2, Double.valueOf(ssActual.getQ2Str()), ssActual.q2Std);
			Assert.assertEquals(ssExpected.rmse, Double.valueOf(ssActual.getRmseStr()), ssActual.rmseStd); 
			Assert.assertEquals(ssExpected.mae, Double.valueOf(ssActual.getMaeStr()), ssActual.maeStd);
		}
	}

	private static void testApplierExport(ModelApplier applier) throws Exception
	{
		String format = "csv";
		File f = null;
		try 
		{
			ExportableSetConfiguration exportConf = new ExportableSetConfiguration();
			exportConf.potentialColumnNames.add("DESCRIPTORS");
			ExportThread thread = applier.getExportThread(exportConf, format);
			Globals.restartAllTransactions(true);
			thread.start();
			thread.join();
			f = new File(thread.eAction.getFullFilePath());
			Assert.assertTrue(f.exists());

			//TODO: Put meaningfull assertions here

		} finally
		{
			if (f != null && f.exists())
				f.delete();
		}
	}

	private static void testApplier(Model model, String testModelFolder) throws Exception
	{
		Basket testSet = OCHEMTestHelpers.generateBasketFromFile(PREFIX, null, getTestsFolder() + "/validation-set.csv");
		Globals.restartAllTransactions(true);

		ModelApplier applier = new ModelApplier();

		applier.compoundsProvider.setBasket(testSet);
		applier.addModel(model);

		applier.defaultTaskPriority = TaskPriority.EXTRA_HIGH;

		applier.start();
		while ( ! applier.isReady()) {
			Thread.sleep(3000);
			applier.update();
			Globals.restartAllTransactions(true);
		}

		testApplierExport(applier);

		ModelApplierTaskProcessor modelTask = applier.modelTasks.get(0);
		modelTask.initialiseModel();
		DataTable dtPredictions = modelTask.wndResult.ports.get(0);
		Assert.assertFalse(modelTask.isError());

		String predictionsPath = getTestsFolder() + "/" + testModelFolder + "/predictions.txt";
		if (generateReferenceFiles)
			savePredictions(predictionsPath, dtPredictions);
		else
			comparePredictions(predictionsPath, dtPredictions);

	}

	/**
	 * Save the predictions to a plain text file, one line per prediction
	 */
	private static void savePredictions(String fileName, DataTable dtPredictions) throws FileNotFoundException
	{
		logger.info("Saving predictions to " + fileName);
		PrintWriter pw = new PrintWriter(new File(fileName));
		dtPredictions.reset();
		while (dtPredictions.nextRow())
			pw.println(dtPredictions.getValue());

		pw.flush();
		pw.close();
	}

	/**
	 * Compare the predictions returned by model applier to the saved predictions
	 */
	private static void comparePredictions(String fileName, DataTable dtPredictions) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(fileName));
		String line;
		dtPredictions.reset();
		int inconsistencies = 0;
		String message = null;
		while ((line = reader.readLine()) != null) {
			dtPredictions.nextRow();

			String expected = NumericalValueStandardizer.getSignificantDigitsStr(Double.parseDouble(line), 2);
			String actual = NumericalValueStandardizer.getSignificantDigitsStr(Double.valueOf("" + dtPredictions.getValue()), 2);
			String difference = NumericalValueStandardizer.getSignificantDigitsStr(Math.abs(Double.parseDouble(expected) - Double.parseDouble(actual)),3);

			if (!expected.equals(actual))
			{
				inconsistencies++;
				message = String.format("Expected: %-4s, Actual: %-4s, Difference: %-4s", expected, actual, difference);
				logger.info(message);
			}
		}
		reader.close();

		if (inconsistencies > INCONSISTENCY_TOLERACE * dtPredictions.getRowsSize())
			throw new AssertionError(String.format("%d inconsistent predictions out of %d. Exemplary error: %s", inconsistencies, dtPredictions.getRowsSize(), message));
	}

	/**
	 * Regenerate all the test case files by rerunning the whole modeling workflow
	 * @throws Exception
	 */
	private static void regenerateTestCases() throws Exception
	{
		ModelCreationTest.generateReferenceFiles = true;

		// Loop over all sub-directories and regenerate the test case files
		String[] childs = new File(getTestsFolder()).list();
		for (String child : childs) {
			if (!child.startsWith(".") && child.contains("oestate") && new File(getTestsFolder() + "/" + child).isDirectory())
				ModelCreationTest.createAndApplyModelTest(child);
		}
	}

	private static String getTestsFolder()
	{
		if (Environment.rootDir == null || Environment.rootDir.equals(""))
			Environment.rootDir = new File(ModelCreationTest.class.getClassLoader().getResource("views.properties").getPath()).getParentFile().getAbsolutePath() + "/../../src/main/webapp";
		return Environment.getWebInfDir() + "/tests/models";
	}

	/**
	 * The main function here is used to generate the new test cases
	 * @param args
	 */
	public static void main(String[] args) 
	{
		new WrapperThread() {

			@Override
			public void wrapped() throws Exception {
				ThreadScope.get().threadClassLoader = ModelStatistics.class.getClassLoader();
				regenerateTestCases();

			}
		}.run();
	}
}

