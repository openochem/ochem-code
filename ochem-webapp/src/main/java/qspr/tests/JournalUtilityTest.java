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

import java.text.ParseException;
import java.text.SimpleDateFormat;

import org.junit.Test;

import qspr.entities.Journal;
import qspr.util.JournalUtility;

@PeriodicTest(schedule = ScheduleType.DAILY, type = "general")
public class JournalUtilityTest
{
	@Test(timeout = 180000)
	public void testSearchJournal() throws Exception
	{
		Journal journal = JournalUtility.fetchJournal("0198-6325"); 
		checkJournalResult(journal);

		journal = JournalUtility.fetchJournal("Medicinal research reviews"); 
		checkJournalResult(journal);

		journal = JournalUtility.fetchJournal("Med Res Rev"); //"1552-8618");
		checkJournalResult(journal);
	}


	@Test(timeout = 180000)
	public void testGetJournalByISSNId() throws Exception
	{
		Journal journal = JournalUtility.fetchJournal("0198-6325");
		journal.update(JournalUtility.fetchJournal("Medicinal research reviews."));
	}

	@Test(timeout = 180000)
	public void testFetchJournalByISSNId() throws Exception
	{
		Journal journal = JournalUtility.fetchJournal("0198-6325"); 
		checkJournalResult(journal);
	}

	@Test(timeout = 180000)
	public void testFetchJournalByTitle() throws Exception
	{
		Journal journal = JournalUtility.fetchJournal("Medicinal research reviews"); 
		checkJournalResult(journal);
	}

	SimpleDateFormat sdf = new SimpleDateFormat("yyyy");

	private void checkJournalResult(Journal journal) throws ParseException
	{
		assertEquals("Medicinal research reviews.", journal.getTitle());
		assertEquals("Med Res Rev", journal.getAbbreviation());
		assertEquals("0198-6325", journal.getISSN());
		assertEquals(sdf.parse("1981"), journal.publish_date);
		assert(journal.publisher.contains("Wiley"));
	}

}
