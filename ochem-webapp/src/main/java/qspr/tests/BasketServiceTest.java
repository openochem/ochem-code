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
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeoutException;

import org.apache.commons.io.FileUtils;
import org.hibernate.criterion.Restrictions;
import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Predicate;
import qspr.entities.Property;
import qspr.entities.PropertyValue;
import qspr.entities.Tag;
import qspr.entities.User;
import qspr.export.ExportableSetConfiguration;
import qspr.util.ExportThread;
import qspr.util.MoleculePeer;

import com.eadmet.business.BasketService;
import com.eadmet.business.DiscretizeBasketOptions;
import com.eadmet.business.DiscretizeBasketRunner;

@PeriodicTest(type = "general", schedule = ScheduleType.DAILY)
public class BasketServiceTest 
{
	@Rule
	public TestRule rule = new LoginRule("BSTest", false);

	private ExperimentalProperty randomRecord(Property p, Article a, ConditionSet conditions, Double value) throws IOException, TimeoutException
	{
		ExperimentalProperty ep = new ExperimentalProperty();
		ep.predicate = Predicate.get("=");
		ep.property = p;
		ep.unit = p.defaultUnit;
		if (value == null)
			ep.value = Math.random();
		else
			ep.value = value;
		ep.article = a;
		ep.molecule = MoleculePeer.getMolecule("CCC");
		ep.owner = ep.introducer = Globals.userSession().user;
		if (conditions != null)
			ep.conditions = conditions.get();
		ep.resolveConflictsAndSave();
		return ep;
	}

	@Test
	@NotifyOnFailure(minConsequentFailures = 3)
	public void basketDiverseTest() throws Exception
	{
		Globals.startAllTransactions();
		Globals.userSession().user.rank = User.SUPER;

		BasketService service = new BasketService();

		//Create some stuff we will need during the test
		String tagName = "Test Tag " + Math.round((Math.random() * 1000));
		Tag tag = new Tag();
		tag.name = tagName;
		tag.type = "property";
		tag.introducer = tag.owner = Globals.userSession().user;
		Globals.session().saveOrUpdate(tag);

		Property p = OCHEMTestHelpers.generateRandomPropertyCondition(Property.TYPE_NUMERIC, false);
		Property c = OCHEMTestHelpers.generateRandomPropertyCondition(Property.TYPE_NUMERIC, true);
		Article a = OCHEMTestHelpers.generateRandomArticle();


		ConditionSet cs = new ConditionSet();
		cs.values.add(new PropertyValue(c, Math.random()));

		List<ExperimentalProperty> eps = new ArrayList<ExperimentalProperty>();

		for (int i=0; i<10; i++)
		{
			Double value = null;

			//For the test with minimal-maximal values in a property distribution
			if (i == 0)
				value = -10D;
			else if (i == 1)
				value = 10D;

			eps.add(randomRecord(p, a, cs, value));
		}

		for (int i=1; i<10; i++)
			eps.get(i).firstEntry = eps.get(0).id; //A bogus primary record linking for the "primary basket" test

		eps.get(0).molecule.mapping1.tags.add(tag);

		//Create our victim basket
		Basket b1 = new Basket();
		b1.user = Globals.userSession().user;
		b1.session = Globals.userSession();
		b1.name = "Test Basket";

		for (ExperimentalProperty ep : eps)
			b1.addEntry(ep);

		Globals.session().saveOrUpdate(b1);
		Globals.restartAllTransactions(true);

		//Test the "edit metadata" that fills the tags, articles and molecules for the basket
		service.fillEditMetadata(b1);
		Assert.assertEquals(1, b1.totalUniqueCompounds.longValue());
		//Assert.assertTrue(b1.tagUsed.contains(tag));
		Assert.assertEquals(a, b1.articleUsed.get(0));

		//Test the list conditions function
		List<Property> conditions = service.listUsedConditions(b1);
		Assert.assertEquals(1, conditions.size());
		Assert.assertEquals(conditions.get(0), c);

		//Test clone and "primary basket". The primary basket requires a better test, but it's complicated
		Basket b2 = service.clone(b1);
		Basket b3 = service.getPrimaryBasket(b1);

		Globals.restartAllTransactions(true);
		b2 = (Basket)Globals.session().get(Basket.class, b2.id);
		b3 = (Basket)Globals.session().get(Basket.class, b3.id);

		for (int i=0; i<b1.entries.size(); i++)
			Assert.assertEquals(b1.entries.get(i).ep.id, b2.entries.get(i).ep.id);

		Assert.assertEquals(1, b3.entries.size());
		Assert.assertEquals(eps.get(0).id, b3.entries.get(0).ep.id); //Our "primary record" for all basket contents

		//Test basket actions  - exclude
		b1 = service.addRecordsFromAnotherBasket(b1, b2, "exclude");
		Globals.restartAllTransactions(true);
		b1 = (Basket)Globals.session().get(Basket.class, b1.id);
		for (int i=0; i<b1.entries.size(); i++)
			Assert.assertTrue(b1.entries.get(i).exclude);

		//Test basket actions  - include
		b1 = service.addRecordsFromAnotherBasket(b1, b2, "include");
		Globals.restartAllTransactions(true);
		b1 = (Basket)Globals.session().get(Basket.class, b1.id);
		for (int i=0; i<b1.entries.size(); i++)
			Assert.assertTrue(!b1.entries.get(i).exclude);

		//Test basket actions  - delete
		b1 = service.addRecordsFromAnotherBasket(b1, b2, "delete");
		Globals.restartAllTransactions(true);
		b1 = (Basket)Globals.session().get(Basket.class, b1.id);
		Assert.assertEquals(0, b1.entries.size());

		//Test basket actions  - add
		b1 = service.addRecordsFromAnotherBasket(b1, b2, "add");
		Globals.restartAllTransactions(true);
		b1 = (Basket)Globals.session().get(Basket.class, b1.id);
		Assert.assertEquals(10, b1.entries.size());


		//Test "fill discretize metadata" - sets properties and their distributions
		service.fillDiscretizeMetadata(b1);
		Assert.assertEquals(p.id, b1.modelProperties.get(0).id);
		Assert.assertEquals(-10D, b1.modelProperties.get(0).distribution.min, 0.01);
		Assert.assertEquals(10D, b1.modelProperties.get(0).distribution.max, 0.01);

		//Test basket export
		String format = "csv";
		File f = null;
		try 
		{
			ExportableSetConfiguration exportConf = new ExportableSetConfiguration();
			exportConf.selectAll = true;
			ExportThread thread = service.getExportThread(b1.id, exportConf, format);
			Globals.restartAllTransactions(true);
			thread.start();
			thread.join();
			f = new File(thread.eAction.getFullFilePath());
			Assert.assertTrue(f.exists());
			String contents = FileUtils.readFileToString(f);
			String[] lines = contents.split("\n");
			Assert.assertEquals(2, lines.length);
			String[] headers = lines[0].split(",");
			Integer propertyIndex = -1;
			for (int i=0; i<headers.length; i++)
				if (headers[i].replaceAll("\"", "").equals(p.getName()+" {measured}"))
					propertyIndex = i;
			Assert.assertTrue(propertyIndex != -1);
			String[] values = lines[1].split(",");
			Assert.assertEquals(-10D, Double.valueOf(values[propertyIndex].replaceAll("\"", "")), 0.001);
		} finally
		{
			if (f != null && f.exists())
				f.delete();
		}

		//Last but not least - test basket discretization. Think of better asserts, add them later
		DiscretizeBasketOptions opts = new DiscretizeBasketOptions();
		opts.basketId = b1.id;
		opts.newBasketName = "Discretized Test Basket";
		opts.newPropertyName = "Discretized Test Property";
		opts.options = new String[]{"test option 1", "test option 2"};
		opts.propertyId = p.id;
		opts.strThresholds = new String[]{"0.5"};
		opts.thresholds = new double[]{0D};

		DiscretizeBasketRunner runner = new DiscretizeBasketRunner(opts);
		runner.start();
		runner.join();

		@SuppressWarnings("unchecked")
		List<Basket> baskets = Globals.session().createCriteria(Basket.class).add(Restrictions.like("name", "Discretized Test Basket")).list();
		Assert.assertEquals(1, baskets.size());
		Basket b = baskets.get(0);
		Assert.assertEquals(2, b.entries.size()); //Our 10 records collapse to 2 unique

		for (BasketEntry be : b.entries)
			Assert.assertEquals("Discretized Test Property", be.ep.property.getName());

		Globals.rollbackAllTransactions();
	}
}
