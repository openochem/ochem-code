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

package qspr.metaserver.util;

import java.io.IOException;
import java.io.Serializable;

import qspr.metaserver.configurations.ModelAbstractConfiguration.ScalingType;
import qspr.workflow.datatypes.DataTable;

public class DataScaling implements Serializable
{

	public double[] slope;
	public double[] bias;

	private static final long serialVersionUID = 1L;

	public DataScaling()
	{
		slope = bias = null;
	}


	public DataScaling(ScalingTable table, ScalingType type, int size) throws IOException
	{

		bias = new double[table.getColumnSize()];
		slope = new double[table.getColumnSize()];

		switch (type)
		{
		case RANGE:
		case RANGE_MINUS1_PLUS1:
			normByRange(table, size, type);
			break;
		case STANDARDIZE:
			standardize(table, size);
			break;
		case NONE:
			break;
		default:
			break;
		}

	}

	public double getInverseSlope(int i)
	{
		return slope[i] == 0 ? 0 : 1. / slope[i];
	}

	public double getInverseBias(int i)
	{
		return slope[i] == 0 ? 0 : -bias[i] / slope[i];
	}

	public double getSlope(int i)
	{
		return slope[i];
	}

	public double getBias(int i)
	{
		return bias[i];
	}

	public void inverseScale(DataTable r)
	{
		if (bias == null)
			return;
		r.reset();
		while (r.hasMoreRows())
		{
			r.nextRow();
			for (int i = 0; i < bias.length; i++)
				r.setValue(i, inverseScale((Double) r.getValue(i), i));
		}
	}

	public double scale(double val, int col)
	{
		if (bias == null)
			return val;
		return bias[col] + slope[col] * val;
	}

	public double inverseScale(double val, int col)
	{
		if (bias == null)
			return val;
		return slope[col] == 0 ? 0 : (val - bias[col]) / slope[col];
	}

	/**
	 * Provide scaling as linear regression
	 * The number of columns and values from the first Table X are used
	 * The number of columns in the second and values in the second table Y should not be less
	 * @param x
	 * @param y
	 */

	public void normByLinearRegression(DataTable x, DataTable y)
	{

		slope = new double[x.getColumnsSize()];
		bias = new double[x.getColumnsSize()];

		for (int i = 0; i < x.getColumnsSize(); i++)
			linearRegression(y, x, x.getRowsSize(), i);

	}

	/**
	 * Performs linear regression and calculates slope and intercept
	 * Y =a+b*X
	 * @param Y
	 * @param X
	 * @param N
	 */
	private void linearRegression(DataTable Y, DataTable X, int N, int col)
	{
		double y[] = new double[N], x[] = new double[N];

		for (int i = 0; i < N; i++)
		{
			y[i] = (Double) Y.getRow(i).getValue(col);
			x[i] = (Double) X.getRow(i).getValue(col);
		}

		double sxy = 0, sx = 0, sy = 0, sxx = 0;

		for (int i = 0; i < N; i++)
		{
			sx += x[i];
			sy += y[i];
			sxy += x[i] * y[i];
			sxx += x[i] * x[i];
		}

		slope[col] = (N * sxy - sx * sy) / (N * sxx - sx * sx);
		bias[col] = (sy - slope[col] * sx) / N;
	}

	private void normByRange(ScalingTable table, int size, ScalingType type) throws IOException
	{
		double min[] = table.getMin(size);
		double max[] = table.getMax(size);

		double maxRange, minRange;

		switch(type) {
		case RANGE:
			maxRange = 1.;
			minRange = 0;
			break;
		case RANGE_MINUS1_PLUS1:
			maxRange = 1;
			minRange = -1;
			break;
		default:
			throw new IOException("This type of normalisation " + type + " is not defined at normByRange");

		}

		for (int i = 0; i < bias.length; i++){
			if(i>=max.length){  // we do not need to normalize above the requested length, i.e. size
				slope[i] = 1;
				bias[0] = 0;
				continue;
			}

			if (max[i] != min[i])
			{
				slope[i] = (maxRange - minRange) / (max[i] - min[i]);
				bias[i] = minRange - min[i] * slope[i];
			}
			else
			{
				slope[i] = 0;
				bias[i] = minRange;
			}
		}

	}

	private void standardize(ScalingTable table, int size) throws IOException
	{
		double mean[] = table.getMean(size);
		double std[] = table.getStd(size);

		for (int i = 0; i < bias.length; i++){

			if(i>=mean.length){  // we do not need to normalize above the requested length, i.e. size
				slope[i] = 1;
				bias[0] = 0;
				continue;
			}

	// incorrect fix 
/*
			if(Math.abs(mean[i]) > std[i] && std[i]/Math.abs(mean[i]) < 10000) {
				slope[i] = 1./mean[i]; // all values will be about one, variation is very small
				bias[i] = 0;
				continue;
			}
*/
			if(std[i]*10<Double.MIN_VALUE)
				std[i] = 1;

			slope[i] = 1. / std[i];
			bias[i] = -mean[i] * slope[i];
		}
	}

}
