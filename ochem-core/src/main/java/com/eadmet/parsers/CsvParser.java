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
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

import au.com.bytecode.opencsv.CSVReader;

import com.eadmet.utils.OCHEMUtils;

public class CsvParser extends SimpleParser
{
	private int importerRow = -1;
	private InputStream globalIs = null;
	private CSVReader globalIm;

	@Override
	public SimpleParser setSource(byte[] data) throws Exception
	{
		super.setSource(data);
		currentRow = -1;
		lastRow = numRecords();
		resolveSchema();
		setCurrentSheet(0);
		return this;
	}

	private int numRecords() throws IOException
	{
		InputStream inp = new BufferedInputStream(new ByteArrayInputStream(data));
		CSVReader reader = new CSVReader(new InputStreamReader(inp));
		try {
			int i=-1;
			String[] values = null;
			while ((values = reader.readNext()) != null)
			{
				if (values.length > 0)
					i++;
				else 
					break;
			}
			return i - 1; //Account for a header which is not a reco
		} finally
		{
			reader.close();
			inp.close();
		}
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
		List<String> currentColumns = new ArrayList<String>();
		columns.add(currentColumns);

		InputStream inp = new BufferedInputStream(new ByteArrayInputStream(data));
		CSVReader reader = new CSVReader(new InputStreamReader(inp));
		try 
		{
			String[] values = reader.readNext();
			for (String value : values)
				currentColumns.add(value);
			return columns;
		} finally
		{
			reader.close();
			inp.close();
		}
	}

	int cachedRowNum = -1;
	List<String> cachedRow = null;

	private List<String> getRowValues(int rowNum) throws IOException
	{
		if (rowNum == cachedRowNum)
			return cachedRow;

		if (importerRow > rowNum)
			createGlobals();

		while (importerRow < rowNum)
		{
			importerRow++;
			globalIm.readNext();
		}

		importerRow++;
		String[] vals = globalIm.readNext();
		List<String> values = new ArrayList<String>();

		if (vals == null)
			return null;

		for (int i = 0; i < vals.length; i++)
			values.add(OCHEMUtils.cleanString(vals[i]));

		boolean emptyRow = true;

		for (String value : values)
			if (value != null)
				emptyRow = false;

		cachedRowNum = rowNum;

		if (emptyRow)
		{
			cachedRow = null;
			return null;
		}
		else
		{
			cachedRow = values;
			return values;
		}
	}

	@Override
	public boolean hasNext() 
	{
		try
		{
			return ((currentRow < lastRow) && (getRowValues(currentRow + 1) != null));
		} catch (IOException e)
		{
			return false;
		}
	}

	@Override
	public List<String> getRow(int rowNumber) throws IOException 
	{
		return getRowValues(rowNumber);
	}

	@Override
	public List<String> next() 
	{
		currentRow++;
		try 
		{
			return getRow(currentRow);
		} catch (Exception e)
		{
			return null;
		}
	}

	private void createGlobals() throws IOException
	{
		importerRow = 0;
		if (globalIm != null)
			globalIm.close();
		if (globalIs != null)
			globalIs.close();

		globalIs = new BufferedInputStream(new ByteArrayInputStream(data));
		int character;
		char special = '\uffff';
		while((character = globalIs.read()) != -1)
			if(character == special)
				throw new IOException("Char " + special + " is reserved and should not be used in the uploaded file");

		globalIs.close();
		globalIs = new BufferedInputStream(new ByteArrayInputStream(data));

		globalIm = new CSVReader(new InputStreamReader(globalIs), au.com.bytecode.opencsv.CSVParser.DEFAULT_SEPARATOR, 
				au.com.bytecode.opencsv.CSVParser.DEFAULT_QUOTE_CHARACTER,  special);
		globalIm.readNext();
	}

	@Override
	public void setCurrentSheet(int currentSheet) throws IOException 
	{
		currentRow = -1;
		currentSheetNumber = currentSheet;
		createGlobals();
		lastRow = numRecords();// - 1;
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
