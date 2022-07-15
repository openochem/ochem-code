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

import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.Tag;

import com.eadmet.business.TagsService;

@PeriodicTest(schedule = ScheduleType.DAILY, type="general")
public class TagsTests
{
	@Rule
	public TestRule loginRule = new LoginRule("TAGS", true);

	@Test
	public void addTest()
	{
		Tag tag = OCHEMTestHelpers.generateRandomTag();
		Basket b = OCHEMTestHelpers.generateRandomBasketNumeric(OCHEMTestHelpers.getRand(), 100);
		@SuppressWarnings("unchecked")
		List<Integer> mp2Ids = Globals.session().createCriteria(BasketEntry.class).add(Restrictions.eq("basket", b)).createAlias("ep", "ep").createAlias("ep.molecule", "mol").createAlias("mol.mapping2", "mp2")
		.setProjection(Projections.groupProperty("mp2.id")).list();

		int uniqueStructureCount = mp2Ids.size();

		TagsService.addTagByMP2(tag, mp2Ids.subList(0, 50));
		Assert.assertEquals(50, tag.getMoleculesCount());

		TagsService.addTagByMP2(tag, mp2Ids.subList(50, mp2Ids.size()));
		Assert.assertEquals(uniqueStructureCount, tag.getMoleculesCount());

	}
}
