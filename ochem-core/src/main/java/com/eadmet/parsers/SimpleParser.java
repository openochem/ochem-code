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
import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.Iterator;
import java.util.List;

import com.eadmet.exceptions.UserFriendlyException;

public abstract class SimpleParser implements Iterable<List<String>>, Iterator<List<String>>
{
	protected byte[] data;
	public String fileName;
	
	public List<String> sheetNames;
	public List<List<String>> sheetColumns;
	
	
	public int currentSheetNumber;
	public int currentRow = 0;
	public int lastRow = 0;
	
	public abstract int getRows();
	
	public abstract boolean hasNext();
	public abstract List<String> next();
	public abstract List<String> getRow(int rowNumber) throws IOException;
	
	public void remove()
	{
		throw new RuntimeException("Unimplemented method");
	}
	
	public abstract void reset();
	
	public abstract List<String> getSheetNames() throws IOException;
	public abstract List<List<String>> getSheetColumns() throws IOException;
	
	public abstract void setCurrentSheet(int sheetNum) throws IOException;
	
	public static SimpleParser getParser(String originalFileName)
	{
		String fileName = originalFileName.toLowerCase();
		SimpleParser parser;	
		if (fileName.endsWith(".xls") || fileName.endsWith(".xlsx"))
			parser = new XlsParser();
		else
		if (fileName.endsWith(".sdf"))
			parser = new SdfParser();
		else
		if (fileName.endsWith(".csv") || fileName.endsWith(".smi"))
			parser = new CsvParser();
		else			
			throw new UserFriendlyException("Provided file has unsupported extension (should be XLS, SDF or CSV)");
		parser.fileName = originalFileName;
		return parser;
	}
	
	public SimpleParser setSource(InputStream is) throws Exception
	{
		BufferedInputStream bis = new BufferedInputStream(is);
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BufferedOutputStream bos = new BufferedOutputStream(baos, 1024);
		byte[] buffer = new byte[1024];
		int bytesRead = -1;
		while ((bytesRead = bis.read(buffer)) != -1)
			bos.write(buffer, 0, bytesRead);
		bos.close();
		bis.close();
		setSource(baos.toByteArray());
		return this;
	}
	
	public SimpleParser setSource(byte[] data) throws Exception
	{
		this.data = data;
		return this;
	}
	
	public SimpleParser setSource(File f) throws Exception
	{
		return setSource(new FileInputStream(f));
	}
	
	protected void resolveSchema() throws IOException
	{
		this.sheetColumns = getSheetColumns();
		this.sheetNames = getSheetNames();
	}
	
	public Iterator<List<String>> iterator()
	{
		reset();
		return this;
	}
	
}
