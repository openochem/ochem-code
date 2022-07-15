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

import java.util.List;

import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.Basket;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Predicate;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.PropertyOptionsFilter;
import qspr.entities.PropertyValue;
import qspr.entities.Tag;
import qspr.entities.Unit;
import qspr.entities.UnitCategory;
import qspr.entities.User;
import qspr.frontend.LabeledValue;
import qspr.util.MoleculePeer;

import com.eadmet.business.PaginationFilter;
import com.eadmet.business.PropertiesAction;
import com.eadmet.business.PropertiesFilter;
import com.eadmet.business.PropertiesFilter.ApprovalStatus;
import com.eadmet.business.PropertiesService;

@PeriodicTest(type = "general", schedule = ScheduleType.DAILY)
public class PropertiesServiceTest 
{
	final String propertyName = "Modified property";
	
	@Rule
	public TestRule rule = new LoginRule("PSTest", false);

	@Test
	@NotifyOnFailure(minConsequentFailures = 3)
	public void propertyActionsTest() throws Exception
	{
		Globals.startAllTransactions();
		Globals.userSession().user.rank = User.SUPER;
		
		PropertiesService service = new PropertiesService();
		
		//Test the stub getters
		List<UnitCategory> uc = service.getUnitCategories();
		List<Predicate> predicates = service.getPredicates();

		
		//Test creating a folder porperty
		PropertiesAction actionFolder = new PropertiesAction();
		actionFolder.name = "Test folder property";
		actionFolder.description = "Test folder property description";
		actionFolder.unitCategoryId = uc.get(0).id;
		actionFolder.defaultUnitId = uc.get(0).getDefaultUnit().id;
		actionFolder.isDirectory = true;
		actionFolder.isCondition = false;
		Property folder = service.edit(actionFolder);

		Assert.assertNotNull(folder.id);
		Assert.assertTrue(folder.isDirectory);

		//Test atomic publish-unpublish operations
		Assert.assertTrue(!folder.isPublished());
		service.publish(folder);
		service.approve(folder);
		Assert.assertTrue(folder.isPublished());
		
		service.unapprove(folder);
		Assert.assertTrue(!folder.approved);
		service.approve(folder);
		Assert.assertTrue(folder.approved);
		service.unapprove(folder);
		Assert.assertTrue(!folder.approved);
		
		
		//Test creating a child property
		PropertiesAction actionChild = new PropertiesAction();
		actionChild.name = "Test" + OCHEMTestHelpers.getRand();
		actionChild.description = "Test property description";
		actionChild.unitCategoryId = uc.get(0).id;
		actionChild.defaultUnitId = uc.get(0).getDefaultUnit().id;
		actionChild.isDirectory = false;
		actionChild.isCondition = false;
		actionChild.isPublic = true;
		actionChild.isApproved = true;
		actionChild.parentId = folder.id;
		Property child = service.edit(actionChild);
		
		//Test removing a child from a folder property
		Assert.assertNotNull(child.id);
		Assert.assertTrue(!child.isDirectory);
		Assert.assertEquals(folder, child.parent);
	
		service.removechild(child);
		Assert.assertNull(child.parent);
		
		service.addchild(folder, child);
		Assert.assertEquals(folder, child.parent);

		service.removechild(child);
		Assert.assertNull(child.parent);
		
		//Test delete property
		service.delete(folder);
		
		//Test create a condition, qualitative
		PropertiesAction actionCondition = new PropertiesAction();
		actionCondition.name = "Test" + OCHEMTestHelpers.getRand();
		actionCondition.description = "Test condition description";
		actionCondition.propertyType = Property.TYPE_QUALITATIVE;
		actionCondition.isDirectory = false;
		actionCondition.isCondition = true;
		Property condition = service.edit(actionCondition);
		Assert.assertTrue(condition.isCondition);
		
		//Test adding removind options
		service.saveOptions(condition, null, null);
		Assert.assertEquals(condition.options.size(), 0);
		
		service.saveOptions(condition, new String[]{"-1","-1","-1"}, new String[]{"TestOption1","TestOption2","TestOption3"});
		Assert.assertEquals(condition.options.size(), 3);
		
		
		//Prepare some needed stuff for further
		Article a = OCHEMTestHelpers.generateRandomArticle();
		
		Basket b = new Basket();
		b.name = "Test basket";
		b.user = Globals.userSession().user;
		Globals.session().saveOrUpdate(b);
		
		
		ConditionSet cs = new ConditionSet();
		PropertyValue pv = new PropertyValue();
		pv.property = condition;
		pv.option = condition.options.get(0);
		cs.values.add(pv);
		
		Unit newUnit = new Unit();
		newUnit.category = uc.get(0);
		newUnit.description = "Temp test unit";
		newUnit.introducer = newUnit.owner = Globals.userSession().user;
		newUnit.setName("TmpUnit");
		Globals.session().saveOrUpdate(newUnit);
		
		//Create some records for property
		for (int i=0; i<5; i++)
		{
			ExperimentalProperty ep = new ExperimentalProperty();
			ep.predicate = predicates.get(0);
			ep.property = child;
			ep.value = i;
			ep.article = a;
			ep.molecule = MoleculePeer.getMolecule("CCC");
			ep.owner = ep.introducer = Globals.userSession().user;
			ep.conditions = cs.get();
			
			if (i < 3)
				ep.unit = child.defaultUnit;
			else
				ep.unit = newUnit;
			
			ep.resolveConflictsAndSave();
			
			b.addEntry(ep);
		}
		
		//Test editing property, adding obligatory conditions
		actionChild.name = propertyName;
		actionChild.obligatoryConditionIds.add(condition.id);
		actionChild.id = child.id;
		child = service.edit(actionChild);
		
		//Test moving property to different unit category
		actionChild.unitCategoryId = uc.get(1).id;
		actionChild.confirmed = true;
		child = service.edit(actionChild);
		
		//Test getting used conditions for property
		service.fillUsedConditions(child);
		Assert.assertEquals(child.conditionsUsed.get(0), condition);
		
		//Test getting used options for our condition
		PropertyOptionsFilter filter = new PropertyOptionsFilter();
		filter.propertyId = condition.id;
		filter.basketId = b.id;
		List<PropertyOption> options = Property.getOptions(filter);
		Assert.assertEquals(options.get(0), condition.options.get(0));
		
		//Test applying a tag
		String tagName = "Test Tag " + Math.round((Math.random() * 1000));
		Tag tag = new Tag();
		tag.name = tagName;
		tag.type = "property";
		tag.introducer = tag.owner = Globals.userSession().user;
		Globals.session().saveOrUpdate(tag);

		service.applyTags(propertyName, tagName);
		Globals.restartAllTransactions(true);

		child = Property.getById(child.id); // we need to get again after the transaction
		tag = Tag.getByID(tag.id);
		
		Assert.assertTrue(child.tags.contains(tag));
		
		//Test filters and selections
		PropertiesFilter f = new PropertiesFilter();
		f.basketId = b.id;
		f.name = "Modified";
		f.tagId = tag.id;
		List<Property> properties = service.get(f, new PaginationFilter(1, 1));
		Assert.assertEquals(child, properties.get(0));
		
		f = new PropertiesFilter();
		f.query = propertyName;
		properties = service.get(f, new PaginationFilter(1, 1));
		Assert.assertEquals(child, properties.get(0));
		
		f = new PropertiesFilter();
		f.id = child.id;
		properties = service.get(f, new PaginationFilter(1, 1));
		Assert.assertEquals(child, properties.get(0));
		
		f = new PropertiesFilter();
		f.approvalStatus = ApprovalStatus.AWAITING_ONLY;
		f.parentId = child.id;
		List<LabeledValue> values = service.getLabels(f, new PaginationFilter(1, 1));
		Assert.assertEquals(0, values.size());
		
		//Hopefully we're done
		Globals.rollbackAllTransactions();
	}
}
