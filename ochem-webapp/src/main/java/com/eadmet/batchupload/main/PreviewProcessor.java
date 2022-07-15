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

package com.eadmet.batchupload.main;

import java.io.IOException;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Mapping2.MMPAIndex;
import qspr.util.AccessChecker;
import qspr.util.UploadContext;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.batchupload.entityschema.EntitiesRemapping;
import com.eadmet.batchupload.main.RecordPreview.PreviewRecordStatus;
import com.eadmet.batchupload.main.RecordPreview.PreviewUploadAction;
import com.eadmet.batchupload.main.RecordPreview.UploadRecordStatus;
import com.eadmet.batchupload.parsers.BatchUploadParser;

public class PreviewProcessor 
{
	BatchUploadParser parser;
	UploadedFileSchema file;
	EntitiesRemapping remapping;
	UploadContext context;
	RecordBuilder builder;
	BatchUploadEventHandler eh;

	RecordValidator validator = new RecordValidator();

	public PreviewProcessor(BatchUploadParser parser, UploadedFileSchema file,  EntitiesRemapping remapping, BatchUploadEventHandler eh, UploadContext context)
	{
		this.parser = parser;
		this.file = file;
		this.remapping = remapping;
		this.context = context;
		this.builder = new RecordBuilder(context);
		this.eh = eh;
	}

	public UploadPreview getPreview() throws IOException
	{
		context.preview = true;
		UploadPreview preview = new UploadPreview();
		RecordStubBuilder rsbuilder = new RecordStubBuilder();
		int totalNum = 0;

		String speed = "";
		long time = Calendar.getInstance().getTimeInMillis();
		int AVERAGE = 100;

		for (int i = 0; i < file.size(); i++) 
		{

			UploadedSheetSchema sheet = file.get(i);

			if (sheet.ignore)
				continue;

			parser.setCurrentSheet(sheet);

			int total = parser.sparser.getRows() + 1;

			for (UploadedRow row : parser)
			{
				List<RecordStub> stubs = rsbuilder.getRowStubs(row);
				for (int recordNum = 0; recordNum < stubs.size(); recordNum++)
				{
					remapping.updateWith(stubs.get(recordNum));
					String recordKey = sheet.getNumber()+"_"+row.number+"_"+recordNum;
					RecordPreview rp = new RecordPreview();
					rp.id = -(++totalNum);
					//					rp.status = PreviewRecordStatus.valid;
					//					rp.action = PreviewUploadAction.save;
					rp.sheet = sheet.getNumber();
					rp.row = row.number;
					rp.record = recordNum;
					preview.add(rp);
					preview.map.put(recordKey, rp);
					fillRecordPreview(rp, stubs.get(recordNum));
					rp.minimize();
				}

				if(row.number != 0 && (row.number % 100) == 0) {
					speed = " (" + ((Calendar.getInstance().getTimeInMillis()- time)/(10*AVERAGE))/100.+ " s/record)";
					time = Calendar.getInstance().getTimeInMillis();
				}

				eh.rowProcessed("Row " + (row.number + 1) + " out of " + total + " processed" + speed);

				if (context.interruptRequested)
					break;
			}

			if (context.interruptRequested)
				break;
		}
		return preview;
	}

	public void uploadFromPreview(UploadPreview preview) throws IOException
	{
		context.preview = false;
		RecordStubBuilder rsbuilder = new RecordStubBuilder();

		String speed = "";
		long time = Calendar.getInstance().getTimeInMillis();
		int AVERAGE = 100;

		for (int i = 0; i < file.size(); i++) 
		{

			UploadedSheetSchema sheet = file.get(i);

			if (sheet.ignore)
				continue;

			parser.setCurrentSheet(sheet);

			Map<Long,Set<Long>> set = new HashMap<Long,Set<Long>>();

			for (UploadedRow row : parser)
			{
				List<RecordStub> stubs = rsbuilder.getRowStubs(row);
				for (int recordNum = 0; recordNum < stubs.size(); recordNum++)
				{
					remapping.updateWith(stubs.get(recordNum));
					String recordKey = sheet.getNumber()+"_"+row.number+"_"+recordNum;
					RecordPreview rp = preview.map.get(recordKey);

					if (rp == null)
						continue;

					if (rp.uploadStatus != null)
						continue;

					saveRecordFromPreview(rp, stubs.get(recordNum),set);
					rp.minimize();
				}

				if(row.number != 0 && (row.number % 100) == 0) {
					speed = " (" + ((Calendar.getInstance().getTimeInMillis()- time)/(10*AVERAGE))/100.+ " s/record)";
					time = Calendar.getInstance().getTimeInMillis();
				}

				eh.rowProcessed("Row " + (row.number + 1) + " out of " + (parser.sparser.lastRow + 1) + " processed" + speed);

				if (context.interruptRequested)
				{
					setUploadStatusSkipped(preview);
					break;
				}
			}

			if (context.interruptRequested)
			{
				setUploadStatusSkipped(preview);
				break;
			}
		}
	}

	private void setUploadStatusSkipped(UploadPreview preview)
	{
		for (RecordPreview rp : preview)
			if (rp.uploadStatus == null)
				rp.uploadStatus = UploadRecordStatus.skipped;
	}

	private void fillRecordPreview(RecordPreview preview, RecordStub stub)
	{
		builder.build(preview, stub);
		validator.validate(preview, context);
		preview.filterMessages();
	}

	private void saveRecordFromPreview(RecordPreview preview, RecordStub stub, Map<Long,Set<Long>> speedSets)
	{
		if (preview.action == PreviewUploadAction.skip)
		{
			preview.uploadStatus = UploadRecordStatus.skipped;
			return;
		}

		fillRecordPreview(preview, stub);

		if (preview.status == PreviewRecordStatus.fatal_error)
		{
			preview.uploadStatus = UploadRecordStatus.fatal_error;
			return;
		}

		try
		{
			switch(preview.action) {
			case merge:
				AccessChecker.requestModificationPermission(preview.duplicate);
				preview.duplicate.merge(preview.ep);

				for (BasketEntry be : preview.ep.basketEntries) 
				{
					Basket b = (Basket)Globals.session().get(Basket.class, be.basket.id);
					BasketEntry beNew = new BasketEntry(preview.duplicate);
					beNew.basket = b;
					b.entries.add(beNew);
				}

				Globals.session().saveOrUpdate(preview.duplicate);
				postprocessSavedRecord(preview.duplicate);
				Globals.session().flush();
				preview.id = preview.duplicate.id;
				preview.uploadStatus = UploadRecordStatus.merged;
				break;

			case put_original_to_basket:
				for (BasketEntry be : preview.ep.basketEntries) 
				{
					Basket b = (Basket)Globals.session().get(Basket.class, be.basket.id);
					be.ep = preview.duplicate;
					if(!speedSets.containsKey(b.id)) {
						Set<Long> set = new HashSet<Long>();
						speedSets.put(b.id, set);
						for(BasketEntry entry:b.entries) 
							set.add(entry.ep.id);
					}

					if (!speedSets.get(b.id).contains(be.ep.id))
					{
						b.entries.add(be);
						speedSets.get(b.id).add(be.ep.id);
						Globals.session().saveOrUpdate(b);
					}
				}
				preview.uploadStatus = UploadRecordStatus.put_original_to_basket;
				preview.ep.basketEntries.clear();
				break;

			case save:
				if (preview.ep.id < 0)
					preview.ep.id = null;

				if (preview.ep.article.id == null)
					Globals.session().saveOrUpdate(preview.ep.article);

				for (BasketEntry be : preview.ep.basketEntries) 
				{
					Basket b = (Basket)Globals.session().get(Basket.class, be.basket.id);
					b.entries.add(be);
				}

				preview.ep.basketEntries.clear();

				if (preview.ep.md5 != null)
					preview.ep.updateHash();

				Globals.session().saveOrUpdate(preview.ep);
				postprocessSavedRecord(preview.ep);
				Globals.session().flush();

				preview.id = preview.ep.id;

				if (preview.status == PreviewRecordStatus.duplicate_external || preview.status == PreviewRecordStatus.duplicate_internal)
					preview.uploadStatus = UploadRecordStatus.saved_duplicate;
				else if (preview.status == PreviewRecordStatus.error)
					preview.uploadStatus = UploadRecordStatus.saved_error;
				else
					preview.uploadStatus = UploadRecordStatus.saved_valid;

				break;
			default:
				break;
			}


		} catch (Exception e)
		{
			preview.uploadStatus = UploadRecordStatus.fatal_error;
			preview.messages.add(BatchUploadMessage.newError(e));
			eh.rowProcessedError("Record exception");
		}
	}

	public void postprocessSavedRecord(ExperimentalProperty ep)
	{
		if (!ep.property.getName().equals(QSPRConstants.DUMMY))
			if (ep.molecule.mapping2.mmpaIndexStatus == MMPAIndex.MMP_IGNORED)
				ep.molecule.mapping2.mmpaIndexStatus = MMPAIndex.MMP_SCHEDULED;
	}

	public void fillRecordPreview(RecordPreview preview) throws IOException
	{
		RecordStubBuilder rsbuilder = new RecordStubBuilder();
		UploadedSheetSchema sheet = file.get(preview.sheet);
		if (!parser.sheet.equals(sheet))
			parser.setCurrentSheet(sheet);
		UploadedRow row = parser.getRow(preview.row);
		List<RecordStub> stubs = rsbuilder.getRowStubs(row);
		RecordStub stub = stubs.get(preview.record);
		remapping.updateWith(stub);
		fillRecordPreview(preview, stub);
	}
}
