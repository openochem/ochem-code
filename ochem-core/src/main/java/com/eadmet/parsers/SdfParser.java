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
import java.util.ArrayList;
import java.util.List;

import qspr.dao.Various;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

public class SdfParser extends SimpleParser
{
	DataTable allData;

	@Override
	public SimpleParser setSource(byte[] data) throws Exception
	{
		InputStream inp = new BufferedInputStream(new ByteArrayInputStream(data));
		try{
			allData = Various.molecule.getAllDataInSDF(inp);
			System.out.println("Source is set");
		}finally{
			inp.close();
		}

		currentRow = -1;
		lastRow = allData.getRowsSize();
		resolveSchema();
		return this;
	}

	@Override
	public List<String> getSheetNames() throws IOException 
	{
		List<String> result = new ArrayList<String>();
		result.add(fileName);
		return result;
	}

	@Override
	public List<List<String>> getSheetColumns() throws IOException 
	{
		List<List<String>> columns = new ArrayList<List<String>>();
		columns.add(allData.getColumns());
		return columns;
	}

	@Override
	public boolean hasNext() 
	{
		return currentRow < lastRow;
	}

	@Override
	public List<String> getRow(int rowNum) 
	{
		if(rowNum >= allData.getRowsSize() || rowNum < 0)
			return null;

		List<String> values = new ArrayList<String>();

		AbstractDataRow row = allData.getRow(rowNum);

		for(int i=0; i< row.size(); i++)
			values.add((String)row.getValue(i));

		return values;	
	}

	@Override
	public List<String> next() 
	{
		currentRow++;
		return getRow(currentRow);
	}

	@Override
	public void setCurrentSheet(int currentSheet) throws IOException 
	{
		currentRow = -1;
		currentSheetNumber = currentSheet;
		lastRow = allData.getRowsSize();
	}

	@Override
	public void reset() 
	{
		currentRow = -1;
	}

	@Override
	public int getRows() {
		return lastRow;
	}
}
