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
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.RuleChain;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import qspr.Globals;
import qspr.export.ExportableSetConfiguration;
import qspr.util.ExportThread;

import com.eadmet.business.DescriptorsStorageService;
import com.eadmet.business.DescriptorsStorageService.CacheUploadResult;
import com.eadmet.descriptorcache.CacheEntry;
import com.eadmet.descriptorcache.DescriptorConfigEntry;
import com.eadmet.descriptorcache.DescriptorsRepository;
import com.eadmet.descriptorcache.DescriptorsRepositoryFactory;

@PeriodicTest(type = "general", schedule = ScheduleType.DAILY)
public class DescriptorsStorageServiceTest 
{
	@Rule
	public TestRule rule = RuleChain.outerRule(new Timeout(50000)).around(new LoginRule("DSSTest", false));

	@Test
	@NotifyOnFailure(minConsequentFailures = 3)
	public void descriptorsCacheManipulateTest() throws Exception
	{
		DescriptorsRepository r = DescriptorsRepositoryFactory.getRepository();
		Globals.startAllTransactions();

		DescriptorsStorageService service = new DescriptorsStorageService();
		// The aggregation of names_size and value_size takes too long.
		// Code is tested for private configs.
		// List<DescriptorConfigEntry> publicConfigs = service.getPublicConfigs();
		// Assert.assertTrue(publicConfigs.size() > 0); 
		List<DescriptorConfigEntry> privateConfigs = service.getPrivateConfigs();
		Assert.assertTrue(privateConfigs.size() == 0);

		DescriptorConfigEntry config = DescriptorsRepositoryTest.createTestConfigEntry("test-deletable-config");
		config.setUser(Globals.userSession().user.login);
		r.saveConfig(config);
		service.deleteCache(config.objectID, Globals.userSession().user.login);
		config = r.getDescriptorConfigById(config.objectID);
		Assert.assertNull(config);

		Globals.rollbackAllTransactions();
	}

	@Test
	@NotifyOnFailure(minConsequentFailures = 3)
	public void descriptorsCacheUploadTest() throws Exception
	{
		DescriptorsRepository r = DescriptorsRepositoryFactory.getRepository();
		try
		{
			Globals.startAllTransactions();
			DescriptorsStorageService service = new DescriptorsStorageService();
			CacheUploadResult result = service.descriptorsCacheUpload(DescriptorsStorageServiceTest.class.getClassLoader().getResourceAsStream("descriptorscache/testCacheUpload.csv"), 
					"testCacheUpload.csv", "UploadedTestDecriptor", null);
			Assert.assertEquals(5, result.count.intValue());
			List<DescriptorConfigEntry> configs = r.getConfigurationsByType("UploadedTestDecriptor");
			Assert.assertEquals(1, configs.size());
			List<CacheEntry> entries = r.getDescriptors(configs.get(0));
			Assert.assertEquals(5, entries.size());
			Globals.rollbackAllTransactions();
		} finally
		{
			r.clearCache("UploadedTestDecriptor", Globals.userSession().user.login);
		}
	}

	@Test
	@NotifyOnFailure(minConsequentFailures = 3)
	public void descriptorsCacheExportTest() throws Exception
	{
		String prefix = "download_test";
		DescriptorConfigEntry config = DescriptorsRepositoryTest.createTestConfigEntry(prefix);
		DescriptorsRepository r = DescriptorsRepositoryFactory.getRepository();
		r.saveConfig(config);
		List<CacheEntry> entries = DescriptorsRepositoryTest.createTestCacheEntries(config, 10, 5, prefix, false);
		r.saveDescriptors(entries);
		String format = "csv";
		File f = null;
		try 
		{
			Globals.startAllTransactions();
			DescriptorsStorageService service = new DescriptorsStorageService();	
			ExportableSetConfiguration exportConf = new ExportableSetConfiguration();
			exportConf.potentialColumnNames.add("DESCRIPTORS");
			ExportThread thread = service.getExportThread(config.objectID, prefix, exportConf, format);
			Globals.restartAllTransactions(true);
			thread.start();
			thread.join();
			f = new File(thread.eAction.getFullFilePath());
			Assert.assertTrue(f.exists());
			String contents = FileUtils.readFileToString(f);
			String[] lines = contents.split("\n");
			Assert.assertEquals(11, lines.length);
			String[] descs = lines[1].split(",");
			Double value = Double.valueOf(descs[0].replace("\"",""));
			float[] values = entries.get(0).getValues();
			Assert.assertEquals(values[0], value, 0.001);

			Globals.rollbackAllTransactions();
		} finally
		{
			r.clearCache(config);
			if (f != null && f.exists())
				f.delete();
		}
	}
}
