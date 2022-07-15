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
import static org.junit.Assert.assertNotNull;

import java.text.SimpleDateFormat;

import org.junit.Test;

import qspr.Globals;
import qspr.entities.Article;
import qspr.util.ISBNUtility;


@PeriodicTest(schedule = ScheduleType.DAILY, type = "general")
public class ISBNUtilityTest
{
	@Test
	public void testGetArticleByISBN() throws Exception
	{
		Globals.startAllTransactions();
		Article article = ISBNUtility.getArticleByISBN("0596009208", new Article());
		assertNotNull(article);
		System.out.println();
		assertEquals("Head first Java", article.getTitle());
		assertEquals("book", article.mediaType);
		assertEquals("Sebastopol, Calif. ; O'Reilly, 2005, c2003.", article.publisher);
		SimpleDateFormat sdf = new SimpleDateFormat("yyyy");
		assertEquals("2005", sdf.format(article.publicationDate));
		Globals.rollbackAllTransactions();
	}
}
