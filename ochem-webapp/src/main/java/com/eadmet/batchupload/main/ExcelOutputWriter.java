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

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;

import qspr.util.UploadContext;

import com.eadmet.batchupload.entityschema.EntitiesRemapping;
import com.eadmet.batchupload.parsers.BatchUploadParser;

public class ExcelOutputWriter 
{
	final static int MAXCOLUMN=256;

	public byte[] generateRemappedExcel(UploadedFileSchema file, EntitiesRemapping remapping, UploadPreview preview, BatchUploadParser parser, BatchUploadEventHandler eh, UploadContext context) throws IOException
	{
		Workbook wb = new HSSFWorkbook();

		Map<String, List<RecordPreview>> rowMap = new HashMap<String, List<RecordPreview>>();
		int cellShift = 0;
		if (preview != null)
		{
			for (RecordPreview p : preview) 
			{
				String key = p.sheet+"_"+p.row;
				if (rowMap.get(key) == null)
					rowMap.put(key, new ArrayList<RecordPreview>());
				rowMap.get(key).add(p);
			}
			cellShift = 3; 
		}

		for (int i = 0; i < file.size(); i++) 
		{
			UploadedSheetSchema sheet = file.get(i);
			parser.setCurrentSheet(sheet);

			Sheet wbsheet = wb.createSheet(sheet.name);
			Row header = wbsheet.createRow(0);

			if (preview != null)
			{
				header.createCell(0).setCellValue("Record status");
				header.createCell(1).setCellValue("Upload status");
				header.createCell(2).setCellValue("Messages");
			}

			int col  = 0;
			for (UploadedColumnSchema column : sheet.getColumns()) {
				if (!column.hidden && (column.realNumber != null))
					header.createCell(column.realNumber + cellShift).setCellValue(column.name);
				if(++col >= MAXCOLUMN) break;
			}

			int rowNum = 0;
			for (UploadedRow row : parser)
			{
				rowNum++;
				Row wbrow = wbsheet.createRow(rowNum);


				if (preview != null && rowMap.get(i+"_"+row.number) != null)
				{
					StringBuilder status = new StringBuilder();
					StringBuilder uploadStatus = new StringBuilder();
					StringBuilder messages = new StringBuilder();

					for (RecordPreview p : rowMap.get(i+"_"+row.number)) 
					{
						if (p.status != null)
							status.append("; "+p.status);
						else
							status.append("; not status");

						if (p.uploadStatus != null)
							uploadStatus.append("; "+p.uploadStatus);
						else
							uploadStatus.append("; not uploaded yet");

						messages.append("; ");
						for (int j = 0; j < p.messages.size(); j++)
						{
							if (j > 0)
								messages.append(", ");
							messages.append(p.messages.get(j).message);
						}
					}
					wbrow.createCell(0).setCellValue(status.substring(2, status.length()));
					wbrow.createCell(1).setCellValue(uploadStatus.substring(2, uploadStatus.length()));
					wbrow.createCell(2).setCellValue(messages.substring(2, messages.length()));
				}

				col  = 0;
				for (UploadedColumnSchema column : sheet.getColumns()) {
					if (!column.hidden && (column.realNumber != null))
					{
						System.out.println(column.realNumber);
						Double numValue = null; String strValue = null;
						try {
							strValue = row.values == null? null : row.values.get(column.realNumber);
							if(strValue != null)
								numValue = Double.valueOf(strValue);
						} catch (Exception e)
						{

						}
						if (numValue != null)
							wbrow.createCell(column.realNumber + cellShift).setCellValue(numValue);
						else
							wbrow.createCell(column.realNumber + cellShift).setCellValue(strValue);
					}
					if(++col >= MAXCOLUMN) break;
				}

				eh.rowProcessed("Row " + (row.number + 1) + " out of " + (parser.sparser.lastRow + 1) + " processed");

				if (context.interruptRequested)
					break;
			}

			if (context.interruptRequested)
				break;

		}

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BufferedOutputStream bos = new BufferedOutputStream(baos);
		wb.write(bos);
		bos.flush();
		bos.close();
		return baos.toByteArray();
	}
}
