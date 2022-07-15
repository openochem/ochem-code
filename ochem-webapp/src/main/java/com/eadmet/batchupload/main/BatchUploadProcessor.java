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

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.business.BatchUpoadBrowserFilter;
import qspr.entities.Property;
import qspr.util.UploadContext;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.batchupload.entityschema.EntitiesRemapping;
import com.eadmet.batchupload.main.UploadedColumnSchema.UploadedColumnType;
import com.eadmet.batchupload.parsers.BatchUploadParser;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.useractions.BatchUpload30Action;
import com.eadmet.useractions.EventFactory;
import com.eadmet.utils.OCHEMUtils;

public class BatchUploadProcessor 
{

	//	public BatchUpload batchUpload;

	public BatchUploadParser parser;

	private UploadedFileSchema file;
	private EntitiesRemapping remapping;
	private UploadPreview preview;


	public UploadContext context = new UploadContext();
	public BatchUploadEventHandler eh = new BatchUploadEventHandler();

	public BatchUploadProcessor setBytes(byte[] data, String fileName) throws Exception
	{
		parser = BatchUploadParser.getParser(data, fileName);
		context.fileName = fileName;
		return this;
	}

	public BatchUploadProcessor setNamedStream(InputStream is, String fileName) throws Exception
	{

		BufferedInputStream bis = new BufferedInputStream(is);

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BufferedOutputStream bos = new BufferedOutputStream(baos, 1024);
		byte[] buffer = new byte[1024];
		int bytesRead = -1;
		while ((bytesRead = bis.read(buffer)) != -1)
			bos.write(buffer, 0, bytesRead);

		bos.close();
		return setBytes(baos.toByteArray(), fileName);		
	}

	public BatchUploadProcessor setFile(File file) throws Exception 
	{
		FileInputStream fis = null;
		try
		{
			fis = new FileInputStream(file);
			return setNamedStream(fis, file.getName());
		} finally
		{
			if (fis != null)
				fis.close();
		}
	}

	void fillSampleValues(UploadedSheetSchema sheet)
	{
		int sampleRows = 10;

		for (UploadedColumnSchema column : sheet.getColumns())
			column.sampleValues.clear();

		for (UploadedRow  row : parser)
		{
			if (sampleRows-- <= 0 || row == null)
				break;

			for (UploadedColumnSchema column : sheet.getColumns())
				if (column.realNumber != null)
				{
					String sv = null;

					if (row.values != null && row.values.size() > column.realNumber)
						sv = StringUtils.abbreviate(row.values.get(column.realNumber), 40);

					if (sv == null) //Nulls in lists are skipped when marshalled to XML
						sv = "";
					column.sampleValues.add(sv);
				}
		}
	}

	boolean isFileSchemaReady()
	{
		if (file == null)
			return false;

		String newMD5 = OCHEMUtils.getMD5(file.toString());

		if (!newMD5.equals(file.md5))
			return false;

		return true;
	}

	boolean isRemappingReady()
	{
		if (!isFileSchemaReady()) //We need the file schema before starting with remapping
			return false;

		// TODO: Check the data MD5 all together, maybe it changed as well?
		if (remapping == null)
			return false;

		if (!remapping.fileSchemaMD5.equals(file.md5)) // New, or someone remapped the file schema
			return false;

		String newMD5 = OCHEMUtils.getMD5(remapping.toString());

		if (!newMD5.equals(remapping.md5))
			return false;

		return true;
	}

	boolean isPreviewReady()
	{
		if (!isRemappingReady())
			return false;

		if (preview == null)
			return false;

		if (!preview.remappingSchemaMD5.equals(remapping.md5))
			return false;

		String newMD5 = OCHEMUtils.getMD5(preview.toString());

		if (!newMD5.equals(preview.md5))
			return false;

		return true;
	}

	boolean isUploadReady()
	{
		if (!isPreviewReady())
			return false;

		for (RecordPreview p : preview.getRecords())
			if (p.uploadStatus == null)
				return false;

		return true;
	}

	void processFileSchema() throws IOException
	{
		if (file == null)
			file = parser.getSchema();

		for (int i = 0; i < file.size(); i++) 
		{
			UploadedSheetSchema sheet = file.get(i);
			parser.setCurrentSheet(sheet);

			if (sheet.addDummyProperty && !sheet.dummyPropertyAdded)
			{
				sheet.addColumn(UploadedColumnType.property, QSPRConstants.DUMMY).setHidden(true);
				sheet.addColumn(UploadedColumnType.value, "0").setHidden(true);
				sheet.dummyPropertyAdded = true;
			} else if (!sheet.addDummyProperty && sheet.dummyPropertyAdded)
			{
				sheet.removeColumn(sheet.getColumns().size()-1);
				sheet.removeColumn(sheet.getColumns().size()-1);
				sheet.dummyPropertyAdded = false;
			}

			ColumnTypeClassifier.classify(sheet);
			fillSampleValues(sheet);
			new SheetSchemaValidator().validate(sheet);
		}	
		file.md5 = OCHEMUtils.getMD5(file.toString());
	}

	void processRemappingSchema() throws IOException
	{
		context.interruptRequested = false;
		boolean needFullRefresh = false;

		if (!isFileSchemaReady())
		{
			needFullRefresh = true;
			processFileSchema();
		}

		if (remapping == null)
		{
			needFullRefresh = true;
			remapping = new EntitiesRemapping();
		}

		remapping.fileSchemaMD5 = file.md5;

		int AVERAGE = 100;
		String speed = "";
		long time = Calendar.getInstance().getTimeInMillis();

		if (needFullRefresh)
		{
			RecordStubBuilder builder = new RecordStubBuilder();
			remapping.records = 0;
			for (int i = 0; i < file.size(); i++) 
			{
				UploadedSheetSchema sheet = file.get(i);

				if (sheet.ignore)
					continue;

				parser.setCurrentSheet(sheet);

				int total = parser.sparser.getRows() + 1;

				for (UploadedRow  row : parser)
				{
					if (row.values == null)
						break;

					if(row.number != 0 && (row.number % 100) == 0) {
						speed = " (" + ((Calendar.getInstance().getTimeInMillis()- time)/(10*AVERAGE))/100.+ " s/record)";
						time = Calendar.getInstance().getTimeInMillis();
					}

					List<RecordStub> stubs = builder.getRowStubs(row);
					for (RecordStub stub : stubs)
						try
					{
							remapping.records++;
							remapping.updateWith(stub);
					}catch(Throwable e) {
						throw new UserFriendlyException("Error is in row: " + (row.number + 1) + " -- try to correct it. See detaild error: " + e.getMessage());
					}

					eh.rowProcessed("Row " + (row.number + 1) + " out of " + total + " processed " + speed);

					if (context.interruptRequested)
						break;
				}

				if (context.interruptRequested)
					break;

			}
		}
		RemappingValidator.validate(remapping, context);
		remapping.md5 = OCHEMUtils.getMD5(remapping.toString());
	}

	void processUploadPreview() throws IOException
	{
		context.interruptRequested = false;

		if (!isRemappingReady())
			processRemappingSchema();

		if (!remapping.valid(false))
			throw new IOException("Remapping schema is not valid, can not proceed with preview");

		eh.rowProcessed("Saving event");
		EventFactory.document("Batch upload 3.0", BatchUpload30Action.previewStart(remapping), null);
		eh.rowProcessed("Saving event");

		context.clearCaches(true);
		preview = new PreviewProcessor(parser, file, remapping, eh, context).getPreview();
		preview.remappingSchemaMD5 = remapping.md5;
		//Validate Preview
		preview.md5 = OCHEMUtils.getMD5(preview.toString());
	}

	void processUpload() throws IOException
	{
		context.interruptRequested = false;

		context.clearCaches(true);

		eh.rowProcessed("Saving event");
		EventFactory.document("Batch upload 3.0", BatchUpload30Action.uploadStart(remapping, preview), null);
		eh.rowProcessed("Saving event");

		new PreviewProcessor(parser, file, remapping, eh, context).uploadFromPreview(preview);

		eh.rowProcessed("Saving event");
		EventFactory.document("Batch upload 3.0", BatchUpload30Action.uploadFinish(remapping, preview), null);
		eh.rowProcessed("Saving event");
	}

	public UploadedFileSchema getFileSchema() throws IOException
	{
		if (!isFileSchemaReady())
			processFileSchema();

		return file;
	}

	public EntitiesRemapping getRemappingSchema() throws IOException
	{
		if (!isRemappingReady())
			processRemappingSchema();

		return remapping;
	}

	public UploadPreview getUploadPreview() throws IOException
	{
		if (!isPreviewReady())
			processUploadPreview();
		return preview;
	}

	public UploadPreview upload() throws IOException
	{
		if (!isUploadReady())
			processUpload();

		return preview;
	}

	public byte[] getUploadReport() throws IOException
	{
		return new ExcelOutputWriter().generateRemappedExcel(file, remapping, preview, parser, eh, context);
	}

	public List<RecordPreview> getRecordPreviews(int startRecord, int totalRecords, BatchUpoadBrowserFilter filter) throws IOException
	{
		context.clearCaches(false);
		context.preview = true;
		PreviewProcessor pbuilder = new PreviewProcessor(parser, file, remapping, eh, context);
		List<RecordPreview> filteredPreview = preview.getRecords(filter);
		List<RecordPreview> originalPreviews = filteredPreview.subList(startRecord, Math.min(startRecord+totalRecords, filteredPreview.size()));
		List<RecordPreview> result = new ArrayList<RecordPreview>();
		for (RecordPreview preview : originalPreviews)
		{
			RecordPreview newPreview = preview.clone();
			pbuilder.fillRecordPreview(newPreview);
			result.add(newPreview);
		}
		return result;
	}

	public List<RecordPreview> getRecordPreviews(int startRecord, int totalRecords) throws IOException
	{
		return getRecordPreviews(startRecord, totalRecords, null);
	}
}

class ColumnTypeClassifier 
{
	@SuppressWarnings("unchecked")
	public static void classify(UploadedColumnSchema column) 
	{
		if (column.name == null)
		{
			column.type = UploadedColumnType.undefined;
			column.ignore = true;
			return;
		}
		Pattern pt = Pattern.compile("(.*)\\{(.*)\\}(.*)",Pattern.DOTALL);
		Matcher m = pt.matcher(column.name);
		if (m.matches())
		{
			column.subscriptName = m.group(2);
			column.name = m.group(1).trim();
		} else
			column.name = column.name.trim();

		column.type = UploadedColumnSchema.UploadedColumnType.fromString(column.name);

		if (column.type.equals(UploadedColumnSchema.UploadedColumnType.undefined))
			column.type = UploadedColumnSchema.UploadedColumnType.fromString(column.name.replaceAll("\\s+", ""));

		if (column.type.equals(UploadedColumnSchema.UploadedColumnType.undefined)) {
			List<Property> list = Globals.session()
					.createCriteria(Property.class)
					.add(Restrictions.eq("shortName", Property.shortName(column.name)))
					.list();
			if (list.size() > 0) {
				Property p = list.get(0);
				if (p.isCondition)
					column.type = UploadedColumnSchema.UploadedColumnType.condition_and_value;
				else
					column.type = UploadedColumnSchema.UploadedColumnType.property_and_value;
			}
		}

		if (column.ignore == null || column.ignore == false)
			column.ignore = column.type.equals(UploadedColumnSchema.UploadedColumnType.undefined);
	}

	public static void classify(UploadedSheetSchema sheet) 
	{
		for (UploadedColumnSchema column : sheet.getColumns())
			classify(column);
	}
}