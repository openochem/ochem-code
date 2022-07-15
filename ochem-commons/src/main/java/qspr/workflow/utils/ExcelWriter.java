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

package qspr.workflow.utils;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import jxl.Cell;
import jxl.NumberCell;
import jxl.Sheet;
import jxl.Workbook;
import jxl.write.Label;
import jxl.write.Number;
import jxl.write.WritableSheet;
import jxl.write.WritableWorkbook;
import jxl.write.WriteException;
import jxl.write.biff.RowsExceededException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

//Facility that helpls to write datatables to XLS
//Make it a public class on necessity
//Midnighter
public class ExcelWriter
{
	private static transient final Logger logger = LogManager.getLogger(ExcelWriter.class);

	WritableWorkbook workbook;
	WritableSheet sheet;
	FileOutputStream fos;
	int currentColumn = 0;

	public ExcelWriter(String path) throws IOException
	{
		new File("", "");
		File f = new File(path);
		f.createNewFile();
		fos = new FileOutputStream(f);
		workbook = Workbook.createWorkbook(fos);
	}

	public ExcelWriter addWnd(WorkflowNodeData wnd) throws Exception
	{
		for (DataTable dt : wnd.ports) 
		{
			String name = ("DT:" + dt.id).replace(":", "");
			addSheet(name);
			addData(dt);
		}

		return this;
	}


	public ExcelWriter addSheet(String name)
	{
		sheet = workbook.createSheet(name, 1);
		currentColumn = 0;

		return this;
	}

	public ExcelWriter addData(DataTable dt) throws RowsExceededException, WriteException
	{
		return addData(dt, 0);
		/*
		int column = 0, n = 1;
		String name = sheet.getName();
		while (column < dt.getColumnsSize())
		{
			if (column != 0)
				addSheet(name+" "+(++n));
			addData(dt, column);
			column += 255;
		}
		return this;
		 */
	}

	private ExcelWriter addData(DataTable dt, int startColumn) throws RowsExceededException, WriteException
	{
		int i = 0, k = 0;
		for (String column : dt.getColumns()) {
			if (k++ < startColumn)
				continue;
			if (i > 255)
				break;
			sheet.addCell(new Label((i++)+currentColumn, 0, column));
		}
		dt.reset();
		while (dt.nextRow())
		{
			i = 0;
			k = 0;
			if (dt.currentRow >= 65534)
			{
				logger.warn("Number of rows in Excel file exceeded the allowed limit. Trimming the content.");
				break;
			}
			for (int ii = 0; ii < dt.getColumnsSize(); ii++) {
				if (k++ < startColumn)
					continue;
				if (i > 255)
					break;
				Object value = dt.getValue(k - 1);
				if (value instanceof Double)
					sheet.addCell(new Number(i + currentColumn, dt.currentRow + 1, (Double)value));
				else if (value instanceof Integer)
					sheet.addCell(new Number(i + currentColumn, dt.currentRow + 1, (Integer)value));
				else if (value instanceof String)
					sheet.addCell(new Label(i + currentColumn, dt.currentRow + 1, (String)value));

				i++;
			}
		}
		currentColumn += dt.getColumnsSize();
		return this;
	}

	public void finalize() throws Exception
	{
		workbook.write();
		workbook.close();
		fos.close();
	}

	// This method should not really be in this class.
	public DataTable getExcelAsDatatable(String path, int sheetNum) throws Exception
	{
		DataTable dt = new DataTable();
		Workbook workbook = Workbook.getWorkbook(new File(path));
		Sheet sheet = workbook.getSheet(sheetNum);
		for (int col = 0; col < sheet.getColumns(); col++)
		{
			String st = sheet.getCell(col, 0).getContents().trim();
			dt.addColumn(st);
		}

		int row = 1;
		boolean hasMoreData = true;
		while (hasMoreData)
		{
			dt.addRow();
			try{
				for (int col = 0; col < sheet.getColumns(); col++)
				{
					Cell cell = sheet.getCell(col, row);
					if (cell instanceof NumberCell)
						dt.setValue(col, ((NumberCell)cell).getValue());
				}
			} catch (ArrayIndexOutOfBoundsException e)
			{
				hasMoreData = false;
				dt.deleteRow(dt.currentRow);
				row--;
			}
			row++;
		}

		logger.info(""+dt.getRowsSize()+" rows created in datatable");

		return dt;
	}
}

