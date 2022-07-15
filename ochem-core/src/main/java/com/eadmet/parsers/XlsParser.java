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

package com.eadmet.parsers;

import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.DateUtil;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;

import com.eadmet.utils.OCHEMUtils;

public class XlsParser extends SimpleParser
{
	private Workbook wb;
	
	@Override
	public SimpleParser setSource(byte[] data) throws Exception
	{
		super.setSource(data);
		InputStream inp = new BufferedInputStream(new ByteArrayInputStream(data));
	    wb = WorkbookFactory.create(inp);
	    setCurrentSheet(0);
	    inp.close();
	    resolveSchema();
	    return this;
	}
	
	private static String getStringValue(Cell cell)
	{
		String value = null;
		if (cell == null)
			return null;
		switch (cell.getCellType())
		{
            case Cell.CELL_TYPE_STRING:
                value = cell.getRichStringCellValue().getString();
                break;
            case Cell.CELL_TYPE_NUMERIC:
                if (DateUtil.isCellDateFormatted(cell))
                   value = cell.getDateCellValue().toString();
                else
                {
                	Double dvalue = cell.getNumericCellValue();
                	if (Math.abs(dvalue - Math.rint(dvalue)) < 1E-9 && Math.abs(dvalue) >= 1)
                		value = new DecimalFormat("#").format(dvalue);
                	else
                		value =  String.valueOf(cell.getNumericCellValue());
                }
                break;
            case Cell.CELL_TYPE_BOOLEAN:
                value = new Boolean(cell.getBooleanCellValue()).toString();
                break;
            case Cell.CELL_TYPE_FORMULA:
            	try
            	{
            		value = new Double(cell.getNumericCellValue()).toString();
            	} catch (Exception e)
            	{
            		value = cell.getCellFormula() + " is error";
 //           		value = cell.getStringCellValue();
            	}
                break;
            default:
            	try
            	{
            		 value = cell.getRichStringCellValue().getString();
            	} catch (Exception e)
            	{
            		e.printStackTrace();
            		value = null;
            	}
                //Or consider empty ?
		}
		return value;
	}
	
	@Override
	public List<String> getSheetNames() throws IOException 
	{
		List<String> names = new ArrayList<String>();
	    int numSheets = wb.getNumberOfSheets();
	    for (int i=0; i<numSheets; i++)
	    {
	    	Sheet sheet = wb.getSheetAt(i);
	    	Row row = sheet.getRow(0);
	    	
	    	if (row == null)
	    		continue;
	    	
	    	names.add(sheet.getSheetName());
	    }
	    return names;
	}
	
	@Override
	public List<List<String>> getSheetColumns() throws IOException 
	{
		List<List<String>> columns = new ArrayList<List<String>>();
	    int numSheets = wb.getNumberOfSheets();
	    for (int i=0; i<numSheets; i++)
	    {
	    	Sheet sheet = wb.getSheetAt(i);
	    	Row row = sheet.getRow(0);
	    	
	    	if (row == null)
	    		continue;
	    	
	    	List<String> currentColumns = new ArrayList<String>();
	    	
	    	for (int j=0; j<row.getLastCellNum(); j++)
	    	{
	    		Cell cell = row.getCell(j);
	    		currentColumns.add(OCHEMUtils.cleanString(getStringValue(cell)));
	    	}
	    	
	    	columns.add(currentColumns);
	    }
	    return columns;
	}
	
	// FIXME: Can add a one-row-cache here to 2 x speedup the parsing stuff
	
	private List<String> getRowValues(int rowNum)
	{
		List<String> values = new ArrayList<String>();
		
		 if (wb.getSheetAt(currentSheetNumber).getRow(rowNum + 1) == null)
			 return null; //Empty row
		
		 List<String> columns = this.sheetColumns.get(currentSheetNumber);
		 for (int i = 0; i < columns.size(); i++)
		{
			Cell cell =  wb.getSheetAt(currentSheetNumber).getRow(rowNum + 1).getCell(i);
			values.add(OCHEMUtils.cleanString(getStringValue(cell)));
		}
		
		boolean emptyRow = true;
		
		for (String value : values)
			if (value != null)
				emptyRow = false;
		
		if (emptyRow)
			return null;
		else
			return values;
	}
	
	@Override
	public boolean hasNext() 
	{
		return ((currentRow < lastRow) && (getRowValues(currentRow + 1) != null));
	}
	
	@Override
	public List<String> getRow(int rowNumber) 
	{
		return getRowValues(rowNumber);
	}

	@Override
	public List<String> next() 
	{
		currentRow++;
		return getRow(currentRow);
	}

	@Override
	public void setCurrentSheet(int sheetNum) throws IOException 
	{
		currentSheetNumber = sheetNum;
		currentRow = -1;
		lastRow = wb.getSheetAt(currentSheetNumber).getLastRowNum();		
	}

	@Override
	public void reset() 
	{
		currentRow = -1;
	}

	@Override
	public int getRows() {
		return lastRow - 1;
	}
}
