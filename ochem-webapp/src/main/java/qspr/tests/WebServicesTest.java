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

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.File;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.RuleChain;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import qspr.Environment;
import qspr.Globals;
import qspr.ThreadScope;
import qspr.entities.Model;
import qspr.entities.Property;
import qspr.entities.Session;
import qspr.workflow.utils.QSPRConstants
;
import qspr.modelling.applier.PropertyPrediction;
import qspr.services.ModelService;
import qspr.services.ModelSummary;

@PeriodicTest(schedule = ScheduleType.DAILY, type = "general")
public class WebServicesTest
{
	public static final String WS_TEST_PREFIX = "WSTest";

	
	@Rule
	public TestRule chain = RuleChain.outerRule(new Timeout(120000))
									 .around(new LoginRule(WS_TEST_PREFIX, false));
	
	@Test
	public void testWsLogin() throws Exception
	{
		ModelService service = new ModelService();
		String guid = service.login(ThreadScope.get().userSession.user.login, ThreadScope.get().userSession.user.getSecret());
		assertNotNull(guid);
	}
	

	@Test
	public void testWsGetModelStatistics() throws Exception
	{
		Session session = Globals.userSession();

		Model mockModel = createMockModel();

		ModelService service = new ModelService();
		PropertyPrediction[] statistics = service.getModelPredictions(session.guid, mockModel.publicId);
		assertTrue(statistics.length > 0);
		assertTrue(statistics[0].getProperty().startsWith("TestProperty_"));
	}
	
	@Test
	public void testWsGetModelSummary() throws Exception
	{
		Session session = Globals.userSession();
		
		Model mockModel = createMockModel();

		ModelService service = new ModelService();
		ModelSummary[] summary = service.getModelSummary(session.guid, mockModel.publicId);
		assertTrue(summary.length == 2);
		assertTrue(summary[0].getPropertyName().startsWith("testproperty"));
		assertTrue(summary[0].getValidationSetName().equals(QSPRConstants.TRAINING));
		assertTrue(summary[1].getPropertyName().startsWith("testproperty"));
		assertTrue(summary[1].getValidationSetName().equals("validation1"));
	}
	

	private Model createMockModel() throws Exception
	{
		Globals.startAllTransactions();
		Model m = OCHEMTestHelpers.trainAModel("TTT", getTestsFolder(), getTestsFolder() + "/ann-cv-oestate/model.ochem.xml", Property.TYPE_NUMERIC, true, true);
		Globals.commitAllTransactions();
		return m;
	}
	
	private static String getTestsFolder()
	{
		if (Environment.rootDir == null || Environment.rootDir.equals(""))
			Environment.rootDir = new File(ModelCreationTest.class.getClassLoader().getResource("views.properties").getPath()).getParentFile().getAbsolutePath() + "/../../src/main/webapp";
		return Environment.getWebInfDir() + "/tests/models";
	}
}

	
	