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

import static org.junit.Assert.assertEquals;

import java.io.InputStream;

import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;

import qspr.Globals;
import qspr.business.toxalert.AlertsFilter;
import qspr.business.toxalert.ScreeningProcessor;
import qspr.entities.Article;
import qspr.entities.Basket;
import qspr.entities.Molecule;
import qspr.entities.Property;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.util.MoleculePeer;

import com.eadmet.business.AlertsService;
import com.eadmet.business.AlertsUploadReport;
import com.eadmet.exceptions.UserFriendlyException;

@PeriodicTest(type = "general", schedule = ScheduleType.BIHOURLY)
public class ToxalertsTest {
	
	@Rule
	public TestRule loginRule = new LoginRule("ToxAlTest", true);
	
	@Test
	@NotifyOnFailure(minConsequentFailures = 5)
	public void toxalertsTest() throws Exception
	{
		CalculationClient client = new CalculationClient("ToxAlertsTest");
		
		if ( ! client.getSupportedTaskTypes().contains(DescriptorsConfiguration.StructuralAlerts))
			throw new UserFriendlyException("StructuralAlerts task type is not supported by metaserver");
		
		ScreeningProcessor processor = new ScreeningProcessor();
		processor.alertsFilter = new AlertsFilter();
		processor.alertsFilter.approvedOnly = false;
		
		Basket basket = new Basket();
		Molecule molecule = MoleculePeer.getMolecule("C=C");
		basket.addMolecule(molecule);
		
		processor.compoundsProvider.setBasket(basket);
		
		Globals.commitAllTransactions();
		processor.start();
		processor.join();
		
		Globals.startAllTransactions();
		
		if (processor.exception != null)
			throw processor.exception;
		
		Assert.assertTrue("expected more than one alert found, but it was " + processor.compoundsByAlerts.index.size(),
				processor.compoundsByAlerts.index.size() >= 1);
	}
	
	@Test
	@NotifyOnFailure(minConsequentFailures = 5)
	public void alertsUploadTest() throws NumberFormatException, Exception
	{
		InputStream fileStream = ToxalertsTest.class.getClassLoader().getResourceAsStream("tests/toxalerts-test-sample.xls");
		Property property = OCHEMTestHelpers.generateRandomPropertyCondition(0, false);
		Article article = OCHEMTestHelpers.generateRandomArticle();
		AlertsUploadReport report = AlertsService.uploadAlerts(fileStream, "toxalerts-test-sample.xls", true, property, article);

		assertEquals(2, report.successes);
		assertEquals(1, report.dublicates);
		
		AlertsFilter filter = new AlertsFilter();
		filter.introducerId = Globals.userSession().user.id;
		assertEquals(2, filter.filterCriteria(null).list().size());
		
		System.out.println(report);
	}
}
