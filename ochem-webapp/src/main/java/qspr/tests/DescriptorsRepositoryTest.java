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

import java.io.IOException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.Random;

import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.RuleChain;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import qspr.metaserver.transport.NoSqlTransport;

import com.eadmet.descriptorcache.CacheEntry;
import com.eadmet.descriptorcache.DescriptorConfigEntry;
import com.eadmet.descriptorcache.DescriptorsRepository;
import com.eadmet.descriptorcache.DescriptorsRepositoryFactory;

@PeriodicTest(type = "general", schedule = ScheduleType.HOURLY)
public class DescriptorsRepositoryTest 
{
	public int TEST_DESCRIPTORS_COUNT = 1000;
	public int TEST_ENTRIES_COUNT = 100;
	public String LITERAL = "repository_test_new";

	@Rule
	public TestRule rule = RuleChain.outerRule(new Timeout(50000)).around(new LoginRule("DSSTest", false));

	private boolean compareEntries(CacheEntry e1, CacheEntry e2)
	{
		String[] names1 = e1.getNames();
		String[] names2 = e2.getNames();
		float[] values1 = e1.getValues();
		float[] values2 = e2.getValues();
		if ((names1.length != names2.length) || (values1.length != values2.length))
			return false;

		for (int i=0; i<names1.length; i++)
		{
			if (!names1[i].equals(names2[i]))
				return false;
			if (values1[i] != values2[i])
				return false;
		}
		return true;
	}

	public static List<CacheEntry> createTestCacheEntries(DescriptorConfigEntry config, int testEntries, int testDescriptors, String perefix, boolean setMP2) throws IOException
	{
		List<CacheEntry> l = new ArrayList<CacheEntry>();

		for (int i=0; i<testEntries; i++)
		{
			Random r = new Random(100500+i);
			CacheEntry c = new CacheEntry();
			c.user = perefix;
			String[] names = new String[testDescriptors];
			float[] values = new float[testDescriptors];
			for (int j = 0; j < testDescriptors; j++)
			{
				names[j] = perefix+"_"+j;
				values[j] = r.nextFloat();
			}
			c.setNamesAndValues(names, values, false);
			c.setMD5(perefix+"_"+i);
			if (setMP2)
				c.mp2 = i;
			c.dateCreated = new Timestamp(Calendar.getInstance().getTimeInMillis());;
			c.error = (i % 2 == 0) ? null : "Error";
			c.attempts = (i % 2 == 0) ? 0 : i;
			c.config = config;
			l.add(c);
		}

		return l;
	}

	public static DescriptorConfigEntry createTestConfigEntry(String prefix)
	{
		DescriptorConfigEntry c = new DescriptorConfigEntry();
		c.description = prefix;
		c.md5 = prefix;
		c.type = prefix;
		c.setUser(prefix);
		return c;
	}

	@Test
	@NotifyOnFailure(minConsequentFailures = 3)
	public void repositoryTest() throws Exception
	{
		long connectionTimeout = System.nanoTime();
		DescriptorsRepository r = DescriptorsRepositoryFactory.getReattemptingRepository();
		connectionTimeout = System.nanoTime() - connectionTimeout;

		try 
		{
			long testTimeout = System.nanoTime();

			DescriptorConfigEntry originalConfig, fetchedConfig;
			List<DescriptorConfigEntry> fetchedConfigList;
			List<CacheEntry> originalEntries, fetchedEntries;
			// Remove prior remains if any
			r.clearCache(LITERAL, LITERAL);
			fetchedConfig = r.getDescriptorConfig(LITERAL);
			Assert.assertNull(fetchedConfig);

			// Checking for saving and returning objectId
			originalConfig = createTestConfigEntry(LITERAL);
			r.saveConfig(originalConfig);
			Assert.assertNotNull(originalConfig.objectID);

			// Checking for fetching by MD5
			fetchedConfig = r.getDescriptorConfig(LITERAL);
			Assert.assertNotNull(fetchedConfig);
			Assert.assertEquals(originalConfig.objectID, fetchedConfig.objectID);

			// Checking for fetching by objectId
			fetchedConfig = r.getDescriptorConfigById(originalConfig.objectID);
			Assert.assertNotNull(fetchedConfig);
			Assert.assertEquals(originalConfig.objectID, fetchedConfig.objectID);

			// Checking for fetching by type
			fetchedConfigList = r.getConfigurationsByType(LITERAL);
			Assert.assertEquals(1, fetchedConfigList.size());
			Assert.assertEquals(originalConfig.objectID, fetchedConfigList.get(0).objectID);

			// Checking for fetching config by user (based on entries) - so far should be empty
			fetchedConfigList = r.getConfigurations(LITERAL);
			Assert.assertEquals(0, fetchedConfigList.size());

			// Checking for saving entries
			originalEntries = createTestCacheEntries(originalConfig, TEST_ENTRIES_COUNT, TEST_DESCRIPTORS_COUNT, LITERAL, true);
			r.saveDescriptors(originalEntries);
			Assert.assertNotNull(originalEntries.get(0).objectID);
			Assert.assertNotNull(originalEntries.get(TEST_ENTRIES_COUNT - 1).objectID);

			// Checking for fetching entries by md5
			fetchedEntries = r.getDescriptors(new String[]{LITERAL+"_0", LITERAL+"_100500", LITERAL+"_3"}, originalConfig);
			Assert.assertTrue(compareEntries(originalEntries.get(0), fetchedEntries.get(0)));
			Assert.assertNull(fetchedEntries.get(1));
			Assert.assertTrue(compareEntries(originalEntries.get(3), fetchedEntries.get(2)));

			//Check for fetching by mp2
			fetchedEntries = r.getDescriptors(new Integer[]{0, 100500, 3}, originalConfig);
			Assert.assertTrue(compareEntries(originalEntries.get(0), fetchedEntries.get(0)));
			Assert.assertNull(fetchedEntries.get(1));
			Assert.assertTrue(compareEntries(originalEntries.get(3), fetchedEntries.get(2)));

			//Check full fetching by config
			fetchedEntries = r.getDescriptors(originalConfig);
			Assert.assertEquals(TEST_ENTRIES_COUNT, fetchedEntries.size());

			// Checking for fetching by user (based on entries) - should be one config and TEST_ENTRIES_COUNT in it
			fetchedConfigList = r.getConfigurations(LITERAL);
			Assert.assertEquals(1, fetchedConfigList.size());
			Assert.assertEquals(TEST_ENTRIES_COUNT, fetchedConfigList.get(0).entriesCount);

			// Checking clear cache by predefined config
			r.clearCache(originalConfig);
			fetchedConfig = r.getDescriptorConfig(LITERAL);
			Assert.assertNull(fetchedConfig);

			testTimeout = System.nanoTime() - testTimeout;
			// Checking that the test can run in 30 seconds. The normal time is less than 5 seconds.
			Assert.assertTrue("Is test run time less than 30 seconds (real:"+testTimeout/1E9+")", testTimeout < 30 * 1E9);

		} catch(AssertionError ee){
			r.clearLostEntries();
			throw ee;
		}
		finally
		{
			r.clearCache(LITERAL, LITERAL);
		}

	}

	public static void main(String[] args) throws Exception
	{
		NoSqlTransport.host = "eco";
		DescriptorsRepositoryTest t = new DescriptorsRepositoryTest();
		t.repositoryTest();
	}
}
