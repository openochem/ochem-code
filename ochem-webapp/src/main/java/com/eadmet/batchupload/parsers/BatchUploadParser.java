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

package com.eadmet.batchupload.parsers;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import com.eadmet.batchupload.main.UploadedColumnSchema;
import com.eadmet.batchupload.main.UploadedFileSchema;
import com.eadmet.batchupload.main.UploadedRow;
import com.eadmet.batchupload.main.UploadedSheetSchema;
import com.eadmet.parsers.SimpleParser;

public class BatchUploadParser implements Iterable<UploadedRow>, Iterator<UploadedRow>
{
	public SimpleParser sparser;
	
	
	public UploadedSheetSchema sheet;
	
	public boolean hasNext()
	{
		return sparser.hasNext();
	}
	
	public UploadedRow next()
	{
		List<String> row = sparser.next();
		return new UploadedRow(sparser.currentRow, row, sheet);
	}
	
	public UploadedRow getRow(int rowNumber) throws IOException
	{
		List<String> row = sparser.getRow(rowNumber);
		return new UploadedRow(rowNumber, row, sheet);		
	}
	
	public void remove()
	{
		throw new RuntimeException("Unimplemented method");
	}
	
	public void reset()
	{
		sparser.reset();
	}
	
	public UploadedFileSchema getSchema() throws IOException
	{
		UploadedFileSchema ufs = new UploadedFileSchema();
		for (int i=0; i<sparser.sheetNames.size(); i++)
		{
			List<String> columns = sparser.sheetColumns.get(i);
			List<UploadedColumnSchema> cs = new ArrayList<UploadedColumnSchema>();
			for (int j = 0; j < columns.size(); j++)
			{
				UploadedColumnSchema ucs = new UploadedColumnSchema(columns.get(j), j);
				cs.add(ucs);
			}
			UploadedSheetSchema uss = new UploadedSheetSchema(sparser.sheetNames.get(i), i, cs);
			ufs.add(uss);
		}
		ufs.setSelectedSheet(0);
		return ufs;
	}
	
	
	public void setCurrentSheet(UploadedSheetSchema sheet) throws IOException
	{
		this.sheet = sheet;
		sparser.setCurrentSheet(sheet.getNumber());
	}

	public static BatchUploadParser getParser(byte[] data, String fileName) throws Exception
	{
		BatchUploadParser parser = new BatchUploadParser();	
		parser.sparser = SimpleParser.getParser(fileName).setSource(data);
		return parser;
	}
	
	protected BatchUploadParser()
	{
	}
	
	public Iterator<UploadedRow> iterator()
	{
		reset();
		return this;
	}
}
