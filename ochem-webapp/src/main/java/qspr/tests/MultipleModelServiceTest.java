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

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.util.List;

import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;

import qspr.Environment;
import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.Property;
import qspr.frontend.MultipleModelsData;

import com.eadmet.business.MultipleModelsService;

@PeriodicTest(type = "general", schedule = ScheduleType.DAILY)
public class MultipleModelServiceTest 
{
	@Rule
	public TestRule rule = new LoginRule("MMSTest", false);

	@Test
	@NotifyOnFailure(minConsequentFailures = 3)
	public void testMultipleModelAggregationAndExport() throws Exception
	{
		final int modelCount = 2;
		Globals.startAllTransactions();
		Property property = OCHEMTestHelpers.generateRandomPropertyCondition(Property.TYPE_NUMERIC, false);
		Basket trainingSet = OCHEMTestHelpers.generateBasketFromFile("MMS_Training", property, getTestsFolder() + "/training-set.csv");
		Basket validationSet = OCHEMTestHelpers.generateBasketFromFile("MMS_Validation", property, getTestsFolder() + "/validation-set.csv");
		for (int i=0; i < 2 * modelCount; i++)
		{
			if (i < modelCount)
				OCHEMTestHelpers.trainAModel("MMS_Model_"+i, trainingSet, validationSet, getTestsFolder() + "/ann-cv-oestate/model.ochem.xml", true, false);
			else
				OCHEMTestHelpers.trainAModel("MMS_Model_"+i, trainingSet, validationSet, getTestsFolder() + "/ann-bagging-oestate/model.ochem.xml", true, false);
		}
		Globals.restartAllTransactions(true);

		MultipleModelsService service = new MultipleModelsService();
		List<MultipleModelsData> data = service.getMultipleModelsData(trainingSet.id, false);
		
		Assert.assertEquals(modelCount, data.get(0).getModelCount()); //Baggings
		Assert.assertEquals(modelCount, data.get(1).getModelCount()); //CVs
		Assert.assertEquals(0, data.get(2).getModelCount()); //No validations
		
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BufferedOutputStream bos = new BufferedOutputStream(baos);
		service.exportAsR(data.toArray(new MultipleModelsData[0]), bos);
		bos.close();
		Assert.assertTrue(baos.size() > 0);

		baos = new ByteArrayOutputStream();
		bos = new BufferedOutputStream(baos);
		service.exportAsXls(data.toArray(new MultipleModelsData[0]), bos);
		bos.close();
		Assert.assertTrue(baos.size() > 0);

		
		Globals.rollbackAllTransactions();
	}
	
	private static String getTestsFolder()
	{
		if (Environment.rootDir == null || Environment.rootDir.equals(""))
			Environment.rootDir = new File(ModelCreationTest.class.getClassLoader().getResource("views.properties").getPath()).getParentFile().getAbsolutePath() + "/../../src/main/webapp";
		return Environment.getWebInfDir() + "/tests/models";
	}
}
