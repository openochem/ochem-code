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

import java.util.ArrayList;
import java.util.List;

import qspr.metaserver.configurations.ConsensusModelConfiguration.ConsensusType;
import qspr.workflow.utils.QSPRConstants;
import qspr.metaserver.util.aggregator.AggregatedValue.SimplePrediction;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

public class AggregatedRowResult
{
	List<AggregatedValue> properties = new ArrayList<AggregatedValue>();
	private String error;
	public boolean keepIndividualValues = false;

	public AggregatedRowResult()
	{

	}

	public void clear() {
		for(AggregatedValue val : properties)
			val.clear();
		error = null;
	}


	public AggregatedRowResult(boolean keepIndividualValues)
	{
		this.keepIndividualValues = keepIndividualValues;
	}

	protected void setRowValues(DataTable dt, int index, int bagIndex, int maxBagIndex) // Add a prediction returned by a calculation server
	{
		AbstractDataRow row = dt.getRow(index);

		if (row.isError())
			error = row.detailedStatus;
		else
			for (int i = 0, k = 0; i < dt.getColumnsSize(); i++)
				if (i == 0 || dt.getColumn(i).toUpperCase().startsWith(QSPRConstants.PREDICTION_RESULT_COLUMN.toUpperCase())) // check if the column is a prediction
				{
					getProperty(k++).addValue((Double)row.getValue(i), bagIndex, maxBagIndex); //add the found value
				}
	}

	public void setRowValues(int property, float[] values, float[] accuracies, ConsensusType type) {

		for (int i = 0; i < values.length; i++)
		{
			if(Float.isNaN(values[i]))
				error = "There is an error with calculation of value: it is NaN";
			else
				getProperty(property).addValue(values[i], accuracies[i], type);
		}		
	}


	protected AbstractDataRow getAggregatedDataRow(ConsensusType type, AbstractDataRow row, boolean allowErrors)
	{
		getProperty(0); // Make sure, there is at least one property for this row

		int n = 0;

		for (AggregatedValue propertyResult : properties)
		{
			if (error != null && !allowErrors)
				row.setError(error);
			else if (propertyResult.count == 0)
				row.setError(error == null ? "This point was not represented in any bag. Bad luck." : error);
			else
			{
				SimplePrediction result = propertyResult.getValues(type);
				row.setValue(n++, result.value);
				row.setValue(n++, result.uncertainity);
				if(Double.isNaN(result.uncertainity))
					row.setError("Sub models failed for this molecule and consensus was not calculated.");

				if (keepIndividualValues)
					row.setValue(n++, propertyResult.getValueArrayIndexed());
			}
		}

		return row;
	}

	protected List<String> getColumns(String dm)
	{
		List<String> columns = new ArrayList<String>();
		if (properties.size() > 1)
			for (int i = 0; i < properties.size(); i++)
			{
				columns.add(QSPRConstants.PREDICTION_RESULT_COLUMN + i);
				columns.add("DM" + i + ":" + dm);
				if (keepIndividualValues)
					columns.add(QSPRConstants.INDIVIDUAL_PREDICTIONS + i);
			}
		else
		{
			columns.add(QSPRConstants.PREDICTION_RESULT_COLUMN);
			columns.add("DM:" + dm);
			if (keepIndividualValues)
				columns.add(QSPRConstants.INDIVIDUAL_PREDICTIONS);
		}

		return columns;
	}

	private AggregatedValue getProperty(int i)
	{
		while (properties.size() <= i)
			properties.add(new AggregatedValue(true));

		return properties.get(i);
	}


}