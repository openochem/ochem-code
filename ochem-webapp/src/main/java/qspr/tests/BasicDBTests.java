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
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.entities.Article;
import qspr.entities.Basket;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Property;

import com.eadmet.exceptions.UserFriendlyException;

@PeriodicTest(schedule = ScheduleType.DAILY, type = "general")
public class BasicDBTests
{
	private static final String PREFIX = "bDBt";
	
	@Rule
	public TestRule loginRule = new LoginRule(PREFIX, true);
	
	@Test
	public void testArticle()
	{
		Article a = OCHEMTestHelpers.generateRandomArticle();
		
		assertEquals(a, Article.getByTitle(a.getTitle()));
		
		assertEquals(a, Article.getById(a.id));
		
		assertEquals(null, a.getUrl());
		a.setLink("link.springer.com/article/10.1007%2Fs10822-011-9440-2");
		assertEquals("http://link.springer.com/article/10.1007%2Fs10822-011-9440-2", a.getUrl());
		
		try
		{
			a.doCheckRights();
		} catch (Exception e)
		{
			assertTrue(e instanceof UserFriendlyException);
		}
		
//		try
//		{
//			assertEquals(a, Article.getArticle(a.pmid.toString()));
//		} catch (Exception e)
//		{
//			assertTrue(e instanceof UserFriendlyException);
//		}
		
		try
		{
			Article.getArticle("99999");
		} catch (Exception e)
		{
			assertTrue(e instanceof UserFriendlyException);
		}
		
		assertEquals(a, Article.getByShortTitle(a.getTitle()));
		
		try
		{
			assertEquals(a, Article.getByISBN("978-0-00-000000-2"));
		} catch (Exception e)
		{
			e.printStackTrace();
		}
		
		assertEquals(a, Article.getByTitle(a.getTitle(), true));
		Article newArticle = Article.getByTitle("StrandKorb", true);
		assertEquals("StrandKorb", newArticle.getTitle());
		
	}
	
	@Test
	public void testProperty()
	{
		Property p = OCHEMTestHelpers.generateRandomPropertyCondition(Property.TYPE_NUMERIC, false);
		
		assertFalse(p.approved);
		p.approve();
		assertTrue(p.approved);
		p.unapprove();
		assertFalse(p.approved);
		
		Globals.setMarshallingOption(MarshallingOption.PROPERTY_RECORD_COUNT);
		System.out.println(p.getPropertyCount());
		
		Globals.setMarshallingOption(MarshallingOption.PROPERTY_UNITS);
		System.out.println(p.getUsedUnit());
		
		Globals.setMarshallingOption(MarshallingOption.PROPERTY_UNITCATEGORY);
		System.out.println(p.getUnitCategory());
		
	}
	
	@Test
	public void testExperimentalProperty()
	{
		Basket basket = OCHEMTestHelpers.generateRandomBasketNumeric(PREFIX + "_testBasket", 2);
		ExperimentalProperty ep = basket.entries.get(0).ep;
		
		// setValueWithPredicate(String newValue) begin
		// setValueWithPredicate with null input
		try
		{
			ep.setValueWithPredicate(null);
		} catch (Exception e)
		{
			assertTrue(e instanceof NullPointerException);
		}
		
		// setValueWithPredicate with empty input
		try
		{
			ep.setValueWithPredicate("");
		} catch (Exception e)
		{
			assertTrue(e instanceof UserFriendlyException);
			assertTrue(e.getMessage().startsWith("No value for property \"TestProperty_"));
		}
		
		// setValueWithPredicate with empty input
		try
		{
			ep.setValueWithPredicate("ochem");
		} catch (Exception e)
		{
			assertTrue(e instanceof UserFriendlyException);
			assertTrue(e.getMessage().startsWith("Invalid value for property \"TestProperty_"));
		}
		// setValueWithPredicate(String newValue) end
		
		// isPublic()
		assertFalse(ep.isPublic());
		
		// getReadonly()
		assertNull(ep.getReadonly());
		
		// merge(ExperimentalProperty ep)
		ExperimentalProperty ep2 = basket.entries.get(1).ep;
		ExperimentalProperty mergedEp = ep.merge(ep2);
		assertEquals(mergedEp, ep);
		
		// cloneForReference()
		assertCloneEqualsOriginal(ep, ep.cloneForReference(), false);
		// fullClone()
		assertCloneEqualsOriginal(ep, ep.fullClone(), true);
		
		// getReferencedProperty()
		assertEquals(ep, ep.getReferencedProperty());
		
		// getStringValue()
		assertEquals(Double.valueOf(ep.value).toString(), ep.getStringValue());
		
		// isDeleted()
		assertFalse(ep.isDeleted());
				
		// TODO Rob 23.07.13: add real assert here
		// setDMValue(String dm, String value)
		ep.setDMValue("newTest", "666");
		
		// getIntroducer()
		assertEquals(Globals.userSession().user, ep.getIntroducer());
		
		// getDuplicate()
		assertEquals(ep.duplicate, ep.getDuplicate());
		
		// approve() begin
		assertFalse(ep.approved);
		try
		{
			ep.approve();
		} catch (Exception e)
		{
			assertTrue(e instanceof UserFriendlyException);
			assertTrue(e.getMessage().startsWith("You are trying to approve a record for a yet unapproved property TestProperty_"));
		}
		ep.property.approve();
//		try
//		{
//			ep.approve();
//		} catch (Exception e)
//		{
//			assertTrue(e instanceof UserFriendlyException);
//			assertTrue(e.getMessage().startsWith("You are not allowed to approve a record for the property TestProperty_"));
//		}
//		ep.property.moderator = Globals.userSession().user;
		try
		{
			ep.approve();
		} catch (Exception e)
		{
			assertTrue(e instanceof UserFriendlyException);
			assertTrue(e.getMessage().startsWith("You are trying to approve a record for a yet unapproved property TestProperty_"));
		}
		assertTrue(ep.approved);
		// approve() end
	}
	
	private void assertCloneEqualsOriginal(ExperimentalProperty original, ExperimentalProperty clone, boolean isFullClone)
	{
		assertEquals(original.property, clone.property);
		assertTrue(original.value == clone.value);
		assertEquals(original.secondValue, clone.secondValue);
		assertEquals(original.predicate, clone.predicate);
		assertEquals(original.unit, clone.unit);
		assertEquals(original.molecule, clone.molecule);
		assertEquals(original.error, clone.error);
		assertEquals(original.conditions, clone.conditions);
		assertEquals(original.option, clone.option);
		if (isFullClone)
		{	
			assertEquals(original.article, clone.article);
			
			assertEquals(original.artLineNum, clone.artLineNum);
			assertEquals(original.artMolId, clone.artMolId);
			assertEquals(original.artPageNum, clone.artPageNum);
			assertEquals(original.artParagraph, clone.artParagraph);
			assertEquals(original.artTableNum, clone.artTableNum);
		}
		assertEquals(original.moleculenames, clone.moleculenames);
		assertEquals(original.rights, clone.rights);
		
	}
	
}
