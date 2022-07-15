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

import java.util.Arrays;

public class Statistics {
	
	private double maximum;
	private double minimum;
	private double arithmeticMean;
	private double geometricMean;
	private double variance;
	private double standardDeviation;
	private double standardError;
	private double median;
	
	private double vals[];
	
	public Statistics(double[] values) {
		
		if (null != values && values.length > 0) {
			this.vals = values;
			this.minimum = this.vals[0];
			this.maximum = this.vals[0];
			this.arithmeticMean = 0;
			this.geometricMean = 1;
			
			double sorted[] = new double[this.vals.length];

			for (int i = 0; i < this.vals.length; i++) {
				
				this.arithmeticMean = this.arithmeticMean + this.vals[i];
				this.geometricMean = this.geometricMean * this.vals[i];
				
				sorted[i] = this.vals[i];
				
				if (this.vals[i] > this.maximum) {
					this.maximum = this.vals[i];
				}
				if (this.vals[i] < this.minimum) {
					this.minimum = this.vals[i];
				}
			}
			
			this.arithmeticMean = this.arithmeticMean / this.vals.length;
			this.geometricMean = Math.pow(this.geometricMean, (1.0 / this.vals.length));
			
			for (int i = 0; i < this.vals.length; i++) {
				this.variance = this.variance + Math.pow((this.vals[i] - this.arithmeticMean), 2);
			}
			
			this.variance = this.variance / this.vals.length;
			this.standardDeviation = Math.sqrt(this.variance);
			this.standardError = Math.sqrt(this.variance / this.vals.length);
			
			Arrays.sort(sorted);
			
			int numb = this.vals.length / 2;
			
			if (ExperimentalDesignHelper.isEven(this.vals.length)) {
				this.median = (sorted[numb] + sorted[numb -1]) / 2;
			}
			else {
				this.median = sorted[numb];
			}
		}
	}

	/**
	 * @return the maximum
	 */
	public double getMaximum() {
		return this.maximum;
	}

	/**
	 * @return the minimum
	 */
	public double getMinimum() {
		return this.minimum;
	}

	/**
	 * @return the arithmeticMean
	 */
	public double getArithmeticMean() {
		return this.arithmeticMean;
	}
	
	/**
	 * @return the geometricMean
	 */
	public double getGeometricMean() {
		return this.geometricMean;
	}

	/**
	 * @return the variance
	 */
	public double getVariance() {
		return this.variance;
	}

	/**
	 * @return the standardDeviation
	 */
	public double getStandardDeviation() {
		return this.standardDeviation;
	}
	
	/**
	 * @return the standardError
	 */
	public double getStandardError() {
		return this.standardError;
	}

	/**
	 * @return the median
	 */
	public double getMedian() {
		return this.median;
	}

	/**
	 * @return the vals
	 */
	public double[] getVals() {
		return this.vals;
	}
	
	public static double getCorrelationCoefficient(double[] var1, double[] var2) {
		
		double coorelationCoefficient = 0;
		Statistics stat1 = new Statistics(var1);
		Statistics stat2 = new Statistics(var2);
		
		for (int i = 0; i < var1.length; i++) {
			
			coorelationCoefficient = coorelationCoefficient + ((var1[i] - stat1.getArithmeticMean()) * (var2[i] - stat2.getArithmeticMean()));
		}
		
		coorelationCoefficient = coorelationCoefficient / (stat1.getStandardDeviation() * stat2.getStandardDeviation() * var1.length);
		
		return coorelationCoefficient;
	}
	
	
}
