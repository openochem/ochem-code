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

package qspr.export;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.IndexedColors;

import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

/**
 * Excel writer for exporting OCHEM data based on Apache POI library 
 * @author midnighter
 *
 */

public class ExcelExportWriter extends ExportWriter
{
	private static transient final Logger logger = LogManager.getLogger(ExcelExportWriter.class);

	HSSFWorkbook workbook;
	LongSheet worksheet;
	int curRow = 1;
	CellStyle csErrorStyle;

	List<LongSheet> sheets = new ArrayList<LongSheet>();

	@SuppressWarnings("rawtypes")
	@Override
	public void writeRow(List row)
	{
		// Create a new worksheet if necessary
		while (currentMolecule.sheet >= sheets.size())
		{
			sheets.add(new LongSheet(data.sheets.get(sheets.size()), row.size()));
			curRow = 1;
		}
		worksheet = sheets.get(currentMolecule.sheet);

		HSSFRow xlsRow = null;
		for (int i = 0; i < row.size(); i++)
		{
			if (i % 256 == 0)
				xlsRow = worksheet.mySheets.get(i / 256).createRow(curRow);
			Object val = row.get(i);
			HSSFCell cell = xlsRow.createCell(i % 256);
			if (val instanceof Double)
				cell.setCellValue((Double) val);
			else
				if (val != null)
					cell.setCellValue("" + val);
			if (currentMolecule.error != null)
				cell.setCellStyle(csErrorStyle);
		}
		curRow++;
	}

	@Override
	public void initialize()
	{
		workbook = new HSSFWorkbook();
		if (data.sheets.isEmpty())
			data.sheets.add(fileName);

		csErrorStyle = workbook.createCellStyle();
		csErrorStyle.setFillForegroundColor(IndexedColors.RED.getIndex());
		csErrorStyle.setFillPattern(CellStyle.SOLID_FOREGROUND);
	}

	@Override
	public void flush() throws IOException
	{
		workbook.write(os);
		super.flush();
	}

	class LongSheet
	{
		List<HSSFSheet> mySheets = new ArrayList<HSSFSheet>();
		String title;

		public LongSheet(String title, int columns)
		{
			this.title = title;
			while (columns > 0)
			{
				addSheet();
				columns -= 256;
			}
		}

		public void addSheet()
		{
			String st = StringUtils.abbreviate(title.replaceAll("/", "_"), 28) + (mySheets.size() > 0 ? "_"+mySheets.size() : "");
			logger.info("Adding sheet "+st);
			HSSFSheet sheet = workbook.createSheet(st);
			HSSFRow row = sheet.createRow(0);
			for (int i = mySheets.size() * 256; i < mySheets.size() * 256 + 256 && i < columns.size(); i++)
				row.createCell(i - mySheets.size() * 256).setCellValue(columns.get(i));
			mySheets.add(sheet);
		}
	}

	@Override
	public void writeSupplementaryData(String key, Object value) 
	{
		HSSFSheet sheet = workbook.createSheet(key);

		if (value instanceof String) // Other types coming soon
		{
			if (((String)value).length() > 32765)
				value = ((String) value).substring(0, 32765);
			sheet.createRow(0).createCell(0).setCellValue((String)value);
		}
		else if (value instanceof Map)
		{
			int i = 0;
			@SuppressWarnings("rawtypes")
			Map map = (Map) value;
			for (Object mapKey : map.keySet())
			{
				HSSFRow row = sheet.createRow(i++);
				row.createCell(0).setCellValue(mapKey.toString());
				row.createCell(1).setCellValue(map.get(mapKey).toString());
			}
		}
		else if (value instanceof DataTable) {
			DataTable tab = (DataTable) value;
			HSSFRow row = sheet.createRow(0);
			for(int i = 0; i<tab.getColumnsSize();i++)
				row.createCell(i).setCellValue(tab.getColumn(i));

			for(int i = 0; i < tab.getRowsSize(); i++) {
				row = sheet.createRow(1+i);
				for(int col = 0; col<tab.getColumnsSize();col++)
					if(tab.getValue(i, col) != null)
					row.createCell(col).setCellValue((Double)tab.getValue(i, col));
			}

		}

	}

	@Override
	public String getFileExtension() {
		return QSPRConstants.EXCEL;
	}

	@Override
	public String getHttpContentType() {
		return "application/xlsx";
	}
}
