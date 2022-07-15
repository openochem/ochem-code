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

import java.text.SimpleDateFormat;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import qspr.Globals;
import qspr.entities.Article;
import qspr.util.PubMedUtility;

@PeriodicTest(schedule = ScheduleType.DAILY, type = "general")
public class PubMedUtilityTest
{

	@Rule
	public TestRule rule = new LoginRule("MNSTest", false);

	@Rule
	public Timeout globalTimeout = new Timeout(180000);

	@Test
	public void getArticleByPubmedId() throws Exception
	{
		Globals.startAllTransactions();

		Article article = PubMedUtility.getArticleByPubmedId(2013l);

		assertEquals("Rapid infusion of sodium bicarbonate and albumin into high-risk premature infants soon after birth: a controlled, prospective trial.", 
				article.getTitle());
		assertEquals("We conducted a controlled, prospective trial to evaluate the effectiveness of rapidly infusing sodium bicarbonate (NaHCO3) and salt-poor albumin into high-risk, premature infants in the first 2 hours of life. Fifty-three infants, randomized into one of four treatment groups, received 8 ml. per kilogram of a solution containing either (A) glucose in water, (B) salt-poor albumin, (C) NaHCO3, or (D) a combination of albumin and NaHCO3. After the initial infusion, the babies received no colloid or alkali solutions until 4 hours of age. We managed them supportively with warmth, appropriate oxygen administration, isotonic fluid infusion, and close monitoring. Among the infants who received alkali, 14 of 26 acquired the respiratory distress syndrome (RDS), 11 died, and four had intracranial hemorrhage. Among babies who received no alkali, RDS occurred in 11 of 27, 5 died, and none had intracranial hemorrhage. These results do not support the common practice of rapidly infusing NaHCO3 into high-risk, premature infants, and they suggest that the early management of such infants needs renewed critical evaluation.", 
				article.articleAbstract);
		for (int i = 0; i < authors.length; i++)
		{
			assertEquals(authors[i], article.authors.get(i).getPrintedName());
		}
		String art = article.journal.getAbbreviation().replaceAll("\\.", "");
		assertEquals("Am J Obstet Gynecol", art);
		SimpleDateFormat sdf = new SimpleDateFormat("d MMM yyyy");
		assertEquals("1 Feb 1976", sdf.format(article.publicationDate));
		assertEquals("124", article.volume);
		assertEquals("3", article.issue);
		assertEquals("263-7", article.pageNumbers);
		assertEquals("2013", article.pmid.toString());
		assertEquals("", article.affiliation);
		assertEquals(null, article.getUrl());
		assertEquals("", article.doi);
		assertEquals(null, article.comment);

		Globals.rollbackAllTransactions();
	}

	private String[] authors = new String[]{"Bland, RD", "Clarke, TL", "Harden, LB"};
}
