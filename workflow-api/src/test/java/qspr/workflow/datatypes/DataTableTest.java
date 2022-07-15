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

package qspr.workflow.datatypes;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

public class DataTableTest 
{
	@Test
	public void columnMergeTest()
	{
		List<DataTable> tables = new ArrayList<DataTable>();
		for (int table=0; table<5; table++)
		{
			DataTable dt = new DataTable(true);
			for (int column=0; column<10; column++)
				dt.addColumn("table_"+table+"_col_"+column);

			for (int row=0; row<100; row++)
			{
				dt.addRow();
				if (row % (table+2) == 0) //Set some amount of rows to errors, different for each generated table
					dt.getCurrentRow().setError("Error for table "+table);
				else
					for (int column=0; column<dt.getColumnsSize(); column++)
						dt.getCurrentRow().setValue(column, table*(row+column)); //Set some value, different for every cell
			}

			tables.add(dt);
		}

		DataTable dt = tables.get(0);
		for (int table=1; table<5; table++)try{
			dt.mergeColumnsWith(tables.get(table)); //Merge all of them down
		}catch(IOException e){}

		Assert.assertEquals(5 * 10, dt.getColumnsSize());
		Assert.assertEquals("table_4_col_5", dt.getColumn(45));
		Assert.assertEquals(4.0 * (50 + 8), dt.getValue(50, 48)); //Table 4, Row 40, Column 8
	}

	@Test
	public void iteratorTest()
	{
		DataTable dt = new DataTable();
		dt.addColumn("VALUE");
		for (int i=0; i<10; i++)
		{
			dt.addRow();
			dt.setValue(i);
		}
		int i = 0;
		for (AbstractDataRow row : dt)
		{
			Assert.assertEquals(i, row.getValue(0));
			i++;
		}
		Assert.assertEquals(10, i);
	}
}
