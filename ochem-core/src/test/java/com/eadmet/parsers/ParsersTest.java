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

import java.io.InputStream;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class ParsersTest 
{
	
	@Test
	public void parserTest() throws Exception
	{
		String[] files = {"parser-test.sdf", "parser-test.xlsx", "parser-test.xls", "parser-test.csv"};
		Integer[] indices = {5, 1, 9, 12, 3, 5, 6};
		for (String file : files) 
		{
			InputStream is = ParsersTest.class.getClassLoader().getResourceAsStream("parserFiles/"+file);
			if (is == null)
				continue;
			SimpleParser sp = SimpleParser.getParser(file).setSource(is);
			sp.reset();
			sp.setCurrentSheet(0);
			List<String> headers = sp.sheetColumns.get(0);
			Assert.assertEquals("Column2", headers.get(2));
			// Continuous access
			int i = 0;
			for (List<String> row : sp)
			{
				if (row == null)
					break;
				Assert.assertEquals(""+(i + 1), row.get(1));
				Assert.assertEquals(""+((i + 1) / 10D), row.get(4));
				i++;
			}
			// Random access
			for (Integer index : indices)
			{
				List<String> row = sp.getRow(index);
				if (row == null)
					continue;
				Assert.assertEquals(""+(index + 1), row.get(1));
				Assert.assertEquals(""+((index + 1) / 10D), row.get(4));
			}
		}
	}
}
