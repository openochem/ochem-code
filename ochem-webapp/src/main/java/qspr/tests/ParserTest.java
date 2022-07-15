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

import java.io.InputStream;

import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;

import qspr.entities.Article;
import qspr.util.EndNoteParser;
import qspr.util.RISParser;

@PeriodicTest(schedule = ScheduleType.DAILY, type = "general")
public class ParserTest
{
	private static final String PREFIX = "parser";
	
	@Rule
	public TestRule loginRule = new LoginRule(PREFIX, true);
	
	@Test
	public void TestRISParser() throws Exception
	{
		InputStream resourceAsStream = ParserTest.class.getClassLoader().getResourceAsStream("tests/parser/jns_jns85_476.ris");
		Article a = new RISParser().parse(resourceAsStream);
		assertRISArticle(a);
	}

	@Test
	public void TestEndNoteParser() throws Exception
	{
		InputStream is = ParserTest.class.getClassLoader().getResourceAsStream("tests/parser/jns_jns85_476.enw");
		Article a = new EndNoteParser().parse(is);
		assertRISArticle(a);
	}
	
	
	private void assertRISArticle(Article a)
	{
		Assert.assertTrue(a.articleAbstract.startsWith("Adult Fisher 344 rats were subjected"));
		Assert.assertTrue(a.authors.size() == 5);
		Assert.assertEquals("4", a.issue);
//		Assert.assertEquals("0022-3085", a.journal.getISSN());
		Assert.assertEquals("Journal of neurosurgery", a.journal.getTitle());
		Assert.assertEquals("476-481", a.pageNumbers);
		Assert.assertEquals("bloodbrainbarrierbreachfollowingcorticalcontusionintherat", 
							Article.shortTitle(a.getTitle()));
		Assert.assertEquals("Blood-brain barrier breach following cortical contusion in the rat", a.getTitle());
		Assert.assertEquals("85", a.volume);
	}
	
}
