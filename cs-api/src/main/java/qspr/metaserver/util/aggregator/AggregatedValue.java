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
import java.util.Collections;
import java.util.List;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.metaserver.configurations.ConsensusModelConfiguration.ConsensusType;

/**
 * An abstract aggregated value with statistical parameters (standard deviation, mean, etc)
 * @author midnighter
 *
 */
public class AggregatedValue
{
	double sum = 0; //arithmetic sum
	double sum2 = 0; // squared sum
	double sumFactors = 0; // sum of all factors
	public int count = 0; // counts number of values averaged
	private double bestAccuracy = Double.MAX_VALUE;
	private double bestValue = 0;

	/**
	 * Do we need to store all values? Its required for the confidence intervals
	 * Allows to save memory if not needed
	 */
	boolean storeAllValues = false;

	SimplePrediction result = new SimplePrediction();

	int maxIndex = -1;
	List<Double> values = new ArrayList<Double>();
	List<Integer> indices = new ArrayList<Integer>();

	/**
	 * Faster operation compared to re-creation the AggregatedValue
	 */

	public void clear() {
		sum = sum2 = sumFactors = bestValue = count = 0;
		bestAccuracy = Double.MAX_VALUE; maxIndex = -1;
		values.clear();;
		indices.clear();
	}

	/**
	 * @param storeAllValues - Do we need to store all the values for advanced calculations (median, confidence interval)?
	 */
	public AggregatedValue(boolean storeAllValues)
	{
		this.storeAllValues = storeAllValues;
	}

	/**
	 * Add a sample value
	 * @param value
	 */
	public void addValue(double value)
	{
		count++;
		sum += value;
		sum2 += value * value;
		if (storeAllValues)
			values.add(value);
	}

	void addValue(double value, int index, int maxIndex)
	{
		addValue(value);
		if (storeAllValues)
			indices.add(index);
		this.maxIndex = maxIndex; 
	}

	void addValue(float value,float accuracy,ConsensusType type) {

		switch(type){
		case AVERAGE:
			addValue(value);
			break;
		case OPTIMAL:
			if(value != Float.NEGATIVE_INFINITY)
				addValue(value);
			break;				
		case WEIGHTED_AVERAGE:
		case RMSE_WEIGHTED:
			addWeightedValue(value, accuracy);
			break;
		case BEST_MODEL:
			addValueWithHighestAccuracy(value, accuracy);
			break;
		}
	}

	SimplePrediction getValues(ConsensusType type) {
		switch(type){
		case AVERAGE:
		case OPTIMAL:
			result.value = getMean();
			result.uncertainity  = getSTD();
			break;
		case BEST_MODEL:
			result.value = getBestValue();
			result.uncertainity = getBestAccuracy();
			break;

		case WEIGHTED_AVERAGE:
			result.value = getMean();
			result.uncertainity  = getSIGMA();
			break;
		case RMSE_WEIGHTED:
			result.value = getMean();
			result.uncertainity  = getSTDWEIGHTED(); // since we do not have individual accuracies and weights are all time the same ...
			break;
		}
		return result;
	}


	/**
	 * Add a sample value with a weighting factor
	 * @param value
	 * @param factor
	 */
	private void addWeightedValue(double value, double factor)
	{
		if(factor<0)throw new UserFriendlyException("Factor should be only positive > 0 "+factor);
		factor *= factor;
		count++;
		sum += value / factor;
		sum2 += value * value / ( factor * factor );
		sumFactors += 1./factor;

		if (storeAllValues)
			values.add(value);
	}

	private void addValueWithHighestAccuracy(double value, double averageRMSE)
	{
		addValue(value);
		if (averageRMSE < bestAccuracy)
		{
			bestAccuracy = averageRMSE;
			bestValue = value;
		}
	}

	/**
	 * Get the normal (or, if available, weighted) mean value
	 * @return
	 */
	public double getMean()
	{
		if (sumFactors == 0)
			return sum / count;
		else
			return sum / sumFactors;
	}

	/**
	 * Get the "radius" of the 95% confidence interval (used to calculate plus-minus values)
	 * @return
	 */
	public double get95ConfidenceInterval()
	{
		if (values != null && values.size() >= 2)
		{
			Collections.sort(values);
			return 1.0 * (values.get((int) (count * 0.975)) - values.get((int) (count * 0.025))) / 2.;
		}

		return 0;
	}

	/**
	 * Get the most accurate value
	 * @return
	 */
	double getBestValue()
	{
		return bestValue;
	}

	/**
	 * Get the most accurate value
	 * @return
	 */
	double getBestAccuracy()
	{
		return bestAccuracy;
	}

	/**
	 * Get sigma value (estimated STD for weighed average)
	 * @return
	 */
	double getSIGMA()
	{
		if(count == 0 || sumFactors == 0) return Double.NaN;
		return Math.sqrt(1./sumFactors);
	}

	/**
	 * Get the standard deviation of the samples
	 * @return
	 */
	public double getSTD()
	{
		if(count == 1) return 0;
		return Math.sqrt((sum2  - sum * sum / count) / (count-1) );
	}

	/**
	 * Get the standard deviation of the samples
	 * @return
	 */
	public double getSTDWEIGHTED()
	{
		if(count == 1) return 0;
		return Math.sqrt((sum2  - sum * sum / count) / (count-1) ) / sumFactors;
	}


	/**
	 * Get the root mean square of the value
	 * @return
	 */
	public double getRMSE()
	{
		return Math.sqrt(sum2 / count);
	}

	float[] getValueArrayIndexed()
	{

		// corresponds to aggregating values in ConsensusModel when the order is preserved

		if(maxIndex == -1 && values.size() >0){
			float[] res = new float[values.size()];

			for (int i=0; i<values.size(); i++)
				res[i] = values.get(i).floatValue();

			return res;
		}

		float[] res = new float[maxIndex];

		for (int i=0; i<indices.size(); i++)
			res[indices.get(i)] = values.get(i).floatValue();

		for (int i=0; i<maxIndex; i++)
			if (!indices.contains(i))
				res[i] = Float.NaN;

		return res;
	}

	class SimplePrediction{
		public double value;
		public double uncertainity;
	}

}
