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

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.Basket;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Property;
import qspr.util.WrapperThread;

import com.eadmet.batchupload.entityschema.EntitiesRemapping;
import com.eadmet.batchupload.main.MultithreadBatchUploadProcessor;
import com.eadmet.batchupload.main.RecordPreview;
import com.eadmet.batchupload.main.RecordPreview.PreviewUploadAction;
import com.eadmet.batchupload.main.RecordPreview.UploadRecordStatus;
import com.eadmet.batchupload.main.TransactionBatchUploadEventHandler;
import com.eadmet.batchupload.main.UploadPreview;
import com.eadmet.batchupload.main.UploadedFileSchema;

@PeriodicTest(type = "general", schedule = ScheduleType.DAILY)
public class BatchUploadsTest 
{
	@Rule
	public TestRule rule = new LoginRule("BUNTest", false);
	
	
	private void buTest() throws Exception
	{
		String fileName = "batch-upload-test-1.xls";
		String basketName = "Basket "+Math.round(Math.random() * 100000);
		Property p1 = OCHEMTestHelpers.generateRandomPropertyCondition(Property.TYPE_NUMERIC, false);
		Property p2 = OCHEMTestHelpers.generateRandomPropertyCondition(Property.TYPE_NUMERIC, false);;
		Article a = OCHEMTestHelpers.generateRandomArticle();
		Globals.restartAllTransactions(true);
		
		MultithreadBatchUploadProcessor buw = new MultithreadBatchUploadProcessor();
		buw.eh = new TransactionBatchUploadEventHandler();
		
		buw.context.hiddenByDefault = true;
		buw.context.allowArticlePubmedSearch = false;
		buw.context.allowMoleculePubchemSearch = false;
		
		buw.setNamedStream(BatchUploadsTest.class.getClassLoader().getResourceAsStream("batchupload/"+fileName), fileName);
		UploadedFileSchema fileSchema =  buw.getFileSchema();
		fileSchema.get(0).ignore = false;
		
		
		fileSchema.get(0).remapColumnToProperty(2, p1.getName()); 
		
		buw.getRemappingSchema();
		buw.waitToFinish();
		Globals.restartAllTransactions(true);
		EntitiesRemapping remappingSchema = buw.getRemappingSchema();
		remappingSchema.properties.get(1).name = p2.getName();
		remappingSchema.articles.get(0).name = "A"+a.id;
		remappingSchema.baskets.get(0).name = basketName;
		
		byte[] data = buw.getUploadReport();
		Assert.assertNotNull(data);
		
		buw.getUploadPreview();
		buw.waitToFinish();
		Globals.restartAllTransactions(true);
		UploadPreview preview = buw.getUploadPreview();
		List<RecordPreview> l = buw.getRecordPreviews(0, 5);
		Assert.assertTrue("Check property 1",l.get(0).ep.property.equals(p1));
		Assert.assertTrue("Check value 1",l.get(0).ep.value == 2.28);
		
		fileSchema.get(0).remapColumnToProperty(2, p2.getName()); 
		
		buw.getRemappingSchema();
		buw.waitToFinish();
		Globals.restartAllTransactions(true);
		remappingSchema.properties.get(1).name = p2.getName();
		remappingSchema.articles.get(0).name = "A"+a.id;
		
		buw.getUploadPreview();
		buw.waitToFinish();
		Globals.restartAllTransactions(true);
		preview = buw.getUploadPreview();
		l = buw.getRecordPreviews(0, 5);
		Assert.assertTrue("Check property 2",l.get(0).ep.property.equals(p2));
		Assert.assertTrue("Check value 2",l.get(0).ep.value == 2.28);
						
		for (RecordPreview recordPreview : preview) 
			recordPreview.action = PreviewUploadAction.skip;
		
		for (int i=0; i<3; i++)
			preview.get(i).action = PreviewUploadAction.save;
		
		Globals.restartAllTransactions(true);
		buw.upload();
		buw.waitToFinish();
		Globals.restartAllTransactions(true);
		preview = buw.upload();
		
		Assert.assertTrue("Check status 1",preview.get(0).uploadStatus == UploadRecordStatus.saved_valid);
		Assert.assertTrue("Check status 2",preview.get(5).uploadStatus == UploadRecordStatus.skipped);
		
		Basket b = Basket.getBasket(basketName, Globals.userSession());
		Assert.assertEquals(3, b.entries.size());
		ExperimentalProperty ep = b.entries.get(0).ep;
		Assert.assertTrue("Check property",ep.property.equals(p2));
		Assert.assertTrue("Check article",ep.article.equals(a));
		Assert.assertTrue("Check value",ep.value == 2.28);
		
		data = buw.getUploadReport();
		Assert.assertNotNull(data);
	}

	@Test
	@NotifyOnFailure(minConsequentFailures = 3)
	public void newBatchUploadTest() throws Exception
	{
		Globals.startAllTransactions();
		buTest();
		Globals.rollbackAllTransactions();
	}
	
	public void buStressTest() throws Exception
	{
		System.out.println("Here we go");
		List<Thread> tlist = new ArrayList<Thread>();
		for (int i=0; i<30; i++)
		{
			WrapperThread t = new WrapperThread()
			{
				@Override
				public void wrapped() throws Exception
				{
					buTest();
				}
			};
			t.userSession = Globals.userSession();
			t.updateSessionTime = false; 
			t.start();
			tlist.add(t);
		}
		
		for (Thread thread : tlist)
			thread.join();
	}
}
