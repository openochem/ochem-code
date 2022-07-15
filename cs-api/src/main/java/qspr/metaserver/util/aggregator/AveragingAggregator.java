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
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

public class AveragingAggregator implements AbstractAggregator
{
	public  String stdDmName = QSPRConstants.BAGGING_STD;
	private List<AggregatedRowResult> rows = new ArrayList<AggregatedRowResult>();

	ConsensusType averagingType = ConsensusType.AVERAGE;

	/**
	 * Indicate whether we should tolerate predictions if some of the aggregated rows are errors
	 * 
	 * If true, we will just aggregate the valid rows.
	 * If false, we will mark the whole result as an error if at least one sub-result fails.
	 */
	boolean allowErrors = true;

	public AveragingAggregator(int size, Boolean keepIndividualPredictions)
	{
		for (int i = 0; i < size; i++)
			rows.add(new AggregatedRowResult(keepIndividualPredictions == null? false:keepIndividualPredictions));
	}


	public AveragingAggregator(int size, ConsensusType averagingType, Boolean allowErrors, Boolean keepIndividualPredictions)
	{
		this.averagingType = averagingType;
		this.allowErrors = allowErrors == null ? false: allowErrors;
		for (int i = 0; i < size; i++)
			rows.add(new AggregatedRowResult(keepIndividualPredictions));
	}

	public static void substituteWithMaxClass(float predictedValues[],float predictedErrors[],int classes) {

		int clas[] = new int[classes];
		for(int i=0;i<predictedValues.length;i++) {
			int cl = Math.round(predictedValues[i]);
			cl = cl<0?0:cl>classes-1?classes-1:cl;
			clas[cl]++;
		}
		int max = 0;
		for(int i=1;i<clas.length;i++)
			if(clas[i]>clas[max])max = i;

		for(int i=0;i<predictedValues.length;i++) {
			predictedValues[i]=max;
			predictedErrors[i]=0;
			}
	
	}

	@Override
	public DataTable getAggregatedResult(boolean compact)
	{
		DataTable dtAggregated = new DataTable(compact); // should be not compact for bagging and consensus to work
		if (rows.size() > 0)
			dtAggregated.addColumns(rows.get(0).getColumns(stdDmName));
		for (AggregatedRowResult aggregatedRowResult : rows)
		{
			AbstractDataRow row = aggregatedRowResult.getAggregatedDataRow(averagingType, dtAggregated.createRow(), allowErrors);
			dtAggregated.addRow(row);
		}

		return dtAggregated;
	}

	@Override
	public void aggregateRow(int rowToStore, DataTable dataTable, int rowInTable, int bagIndex, int maxBagIndex)
	{
		AggregatedRowResult pr = rows.get(rowToStore); // row to make aggregation
		pr.setRowValues(dataTable, rowInTable, bagIndex, maxBagIndex);
	}

	public void addPredictions(int row,float values[], float accuracies[]) {
		rows.get(row).setRowValues(0, values, accuracies, averagingType);
	}

	public void addPredictions(int row, int property, float values[], float accuracies[]) {
		rows.get(row).setRowValues(property, values, accuracies, averagingType);
	}

	/**
	 * Clean content of all rows 
	 */

	public void cleanRows() {
		for (AggregatedRowResult row: rows)
			row.clear();
	}

}
