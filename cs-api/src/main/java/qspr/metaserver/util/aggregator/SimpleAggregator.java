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

package qspr.metaserver.util.aggregator;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

public class SimpleAggregator implements AbstractAggregator
{
	private static final Logger logger = LogManager.getLogger(SimpleAggregator.class);
	
	int dataSize;
	DataTable dt = null;
	Map<Integer, AbstractDataRow> results = new HashMap<Integer, AbstractDataRow>();

	public SimpleAggregator(int size)
	{
		dataSize=size;
	}

	@Override
	public void aggregateRow(int rowToStore, DataTable dataTable, int rowInTable, int bagIndex, int maxBagIndex)
			throws IOException 
	{
		if (dt == null)
			dt = dataTable.getEmptyCopy();

		if(!dataTable.isCompactRowFormat())
			throw new IOException("The dataTable is non-compact:"
					+ dataTable.getRow(rowInTable));

		results.put(rowToStore, dataTable.getRow(rowInTable));
	}

	@Override
	public DataTable getAggregatedResult(boolean compact)
			throws IOException {

		logger.info("Aggregating " + dataSize + " rows");

		if (dataSize != results.size())
			throw new IOException(
					"There are different sizes of data table in aggregated results:"
							+ results.size() + " and original table:"
							+ dataSize);

		for (int i = 0; i < dataSize; i++)
			dt.addRow(results.get(i));

		return dt;
	}

}
