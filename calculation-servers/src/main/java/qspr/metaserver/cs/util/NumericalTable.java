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

package qspr.metaserver.cs.util;

import java.io.IOException;

import qspr.metaserver.configurations.ModelAbstractConfiguration.ScalingType;
import qspr.metaserver.util.DataScaling;
import qspr.metaserver.util.ScalingTable;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.utils.NumericalValueStandardizer;

/**
 * An abstract wrapper around numerical tables Supports scaling and basic statistics
 * 
 */
abstract public class NumericalTable implements ScalingTable {

	public static final String MISSED_VALUES = "999999";
	public static final double MISSED_VALUE = Double.parseDouble(MISSED_VALUES);
	
	DataTable data;
	DataScaling scale;

	/**
	 * This is the main and recommended function to get values and descriptors.
	 *  It returns String, i.e. correctly formatted values ready to be saved to file descriptors and variables.
	 * 
	 * @param mol
	 * @return
	 * @throws IOException
	 */

	public abstract String[] getScaledValuesString(int molecule) throws IOException;

	final protected double[] getScaledValues(int molecule) throws IOException { // only one value is present for one row in Labels table
		String[] values = getScaledValuesString(molecule);
		double val[] = new double[values.length];
		for(int i =0; i < values.length; i++)
			val[i] = Double.parseDouble(values[i]);
		return val;
	}

	/**
	 * Provides non scaled values for various normalization purposes
	 */


	protected double[] getNonScaledValues(int molecule) throws IOException {
		DataScaling storedScale = scale;
		scale = null;
		double val[]  =  getScaledValues(molecule);
		scale = storedScale;
		return val;
	}


	public void createScaling(ScalingType type, int size) throws IOException {
		scale = new DataScaling(this, type, size);
	}

	@Override
	public int getColumnSize()
	{
		return data.getColumnsSize();
	}

	public int getDataSize() {
		return data == null ? 0 : data.getRowsSize();
	}

	// TODO eliminate code which uses this function!
	public DataTable getRawData() {
		return data;
	}

	public AbstractDataRow getRawRow(int index) {
		return data.getRow(index);
	}

	public DataScaling getScaling() {
		return scale;
	}

	public void setScaling(DataScaling scale) {
		this.scale = scale;
	}

	String scaleValue(double val, int col) {
		return NumericalValueStandardizer.getSignificantDigits(scale == null ? val : scale.scale(val, col));
	}

	protected float floatScaleValue(float val, int col) {
		return (scale == null ? val : (float)scale.scale(val, col));
	}

	/**
	 * Returns std of descriptors and values (including options and qualitative descriptors)
	 * 
	 * @throws IOException
	 */

	public double[] getStd(int size) throws IOException {
		double[] means = getMean(size);
		double std[] = new double[means.length];

		for (int mol = 0; mol < size; mol++) {
			double v[] = getNonScaledValues(mol);
			for (int i = 0; i < std.length; i++)
				std[i] += (v[i] - means[i]) * (v[i] - means[i]);
		}

		for (int i = 0; i < std.length; i++)
			std[i] = (float)Math.sqrt(std[i] / size);

		return std;
	}

	/**
	 * Return means values of descriptors & value
	 */

	public double[] getMean(int size) throws IOException {

		double val[] = getNonScaledValues(0);

		for (int mol = 1; mol < size; mol++) {
			double v[] = getNonScaledValues(mol);
			for (int i = 0; i < val.length; i++)
				val[i] += v[i];
		}
		for (int i = 0; i < val.length; i++)
			val[i] /= size;

		return val;
	}

	public double[] getMin() throws IOException {
		return getMin(getDataSize());
	}

	public double[] getMax() throws IOException {
		return getMax(getDataSize());
	}


	public double[] getMin(int size) throws IOException {

		double val[] = getNonScaledValues(0);

		for (int mol = 0; mol < size; mol++) {
			double v[] = getNonScaledValues(mol);
			for (int i = 0; i < val.length; i++) {
				if (val[i] > v[i])
					val[i] = v[i];
			}
		}
		return val;
	}

	public double[] getMax(int size) throws IOException {

		double val[] = getNonScaledValues(0);

		for (int mol = 0; mol < size; mol++) {
			double v[] = getNonScaledValues(mol);
			for (int i = 0; i < val.length; i++)
				if (val[i] < v[i])
					val[i] = v[i];
		}
		return val;
	}

	public NumericalTable getSlice(int fromIndex, int toIndex) {
		NumericalTable table = getCopy();
		table.data = data.getSlice(fromIndex, toIndex);
		return table;
	}

	public NumericalTable getCopy() {
		try {
			NumericalTable table = this.getClass().newInstance();
			table.data = data;
			return table;
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	@Override
	public String toString() {
		return data == null ? "" : data.toString();
	}
	
	public NumericalTable deleteFailed(boolean skip[]) {
		NumericalTable table = getCopy();
		table.data = data.deleteFailed(skip);
		return table;
	}
	

}
