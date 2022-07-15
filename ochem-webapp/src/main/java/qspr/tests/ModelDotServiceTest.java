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

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.freehep.util.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;

import qspr.Environment;
import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.Property;
import qspr.modelling.ModelStatistics;

import com.eadmet.business.ModelDotService;
import com.eadmet.business.PaginationFilter;

@PeriodicTest(type = "general", schedule = ScheduleType.DAILY)
public class ModelDotServiceTest 
{
	@Rule
	public TestRule rule = new LoginRule("MDSTest", false);

	private static String getTestsFolder()
	{
		if (Environment.rootDir == null || Environment.rootDir.equals(""))
			Environment.rootDir = new File(ModelCreationTest.class.getClassLoader().getResource("views.properties").getPath()).getParentFile().getAbsolutePath() + "/../../src/main/webapp";
		return Environment.getWebInfDir() + "/tests/models";
	}

	@SuppressWarnings("deprecation")
	@Test
	@NotifyOnFailure(minConsequentFailures = 3)
	public void testModelDots() throws Exception
	{
		Globals.startAllTransactions();
		// Prepare models
		Property property = OCHEMTestHelpers.generateRandomPropertyCondition(Property.TYPE_NUMERIC, false);
		Basket trainingSet = OCHEMTestHelpers.generateBasketFromFile("MMS_Training", property, getTestsFolder() + "/training-set.csv");
		Basket validationSet = OCHEMTestHelpers.generateBasketFromFile("MMS_Validation", property, getTestsFolder() + "/validation-set.csv");
		List<Model> models = new ArrayList<Model>();
		for (int i=0; i < 2; i++)
		{
			if (i < 1)
				models.add(OCHEMTestHelpers.trainAModel("MMS_Model_"+i, trainingSet, validationSet, getTestsFolder() + "/ann-cv-oestate/model.ochem.xml", true, false));
			else
				models.add(OCHEMTestHelpers.trainAModel("MMS_Model_"+i, trainingSet, validationSet, getTestsFolder() + "/ann-bagging-oestate/model.ochem.xml", true, false));
		}
		Globals.restartAllTransactions(true);
		//
		ModelDotService service = new ModelDotService();

		for (Model model : models) 
		{
			model = (Model)Globals.session().get(Model.class, model.id); //Re-fetch the model
			// Test model dots - all valid
			ModelMapping mapping = model.modelMappings.get(0);
			ModelStatistics stat = (ModelStatistics)mapping.statisticsOriginal.getObject();

			Long modelMappingId = mapping.id;
			Long statNum = 0L;
			Long epId = stat.sets.get(0).points.get(0).id;

			List<ExperimentalProperty> dots = service.getModelDotList(model, statNum, epId, modelMappingId);

			Assert.assertTrue(dots.size() > 0);
			Assert.assertNotNull(dots.get(0).predicted);

			// Test model dots errors with pagination
			PaginationFilter pager = new PaginationFilter(1, 0);
			List<ExperimentalProperty> errors = service.getErrorList(model, pager, false, null);
			List<Long> ids = ModelDotService.getErrorRecordsIds(model.modelMappings.get(0), false, null);
			Assert.assertEquals(Integer.valueOf(errors.size()).longValue(), pager.totalSize.longValue());
			Assert.assertEquals(errors.size(), ids.size());

			//Test descriptors
			List<String> descriptors = service.getDescriptorList(model, modelMappingId, epId);
			Assert.assertNotNull(descriptors); //Meaningless check... however without a previously prepared model we cannot assert anything meaningfull
		}

		Globals.rollbackAllTransactions();
	}
}
