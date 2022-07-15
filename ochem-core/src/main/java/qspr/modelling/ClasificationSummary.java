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

package qspr.modelling;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.commons.math3.random.MersenneTwister;

import com.eadmet.utils.NumericalValueStandardizer;

@XmlRootElement
public class ClasificationSummary implements Serializable
{
	private static final long serialVersionUID = 1L;

	// TP,TN,FP,FN and MCC are not determined and thus all are 0s for multiclass
	// classification
	// For 2-class classification it is supposed that
	// 0 corresponds to positive class
	// 1 corresponds to negative class

	public Integer tp;
	public Integer tn;
	public Integer fp;
	public Integer fn;

	public RangedValue accuracyTotal = new RangedValue();
	public RangedValue accuracyBalanced = new RangedValue();
	public RangedValue mcc = new RangedValue();
	public RangedValue auc = new RangedValue();

	@XmlElement
	private List<RangedValue> valueByClass = new ArrayList<RangedValue>();

	@XmlElement
	private List<RangedValue> sensitivityByClass = new ArrayList<RangedValue>();

	@XmlTransient
	private Map<String, Node> nodesMap = new HashMap<String, Node>();

	@XmlElement
	private List<Node> nodes = new ArrayList<Node>();

	private Node getNode(Integer real, Integer predicted)
	{
		Node node = nodesMap.get("" + real + "," + predicted);

		if (node == null)
		{
			node = new Node();
			node.real = real * 1L;
			node.predicted = predicted * 1L;
			nodesMap.put("" + real + "," + predicted, node);
			nodes.add(node);
		}

		return node;
	}

	public int getNumberOfClasses() {
		int numOfClasses = 0;
		for (Node node : nodes)
		{
			numOfClasses = Math.max(Math.round(node.real) + 1, numOfClasses);
		}

		return numOfClasses;
	}

	public void calculateStatistics(int replicas)
	{
		int numOfPoints = 0;
		int numOfClasses = 0;
		for (Node node : nodes)
		{
			numOfPoints += node.count.intValue();
			numOfClasses = Math.max(Math.round(node.real) + 1, numOfClasses);
		}

		int[] real = new int[numOfPoints + 1], predicted = new int[numOfPoints + 1];

		int k = 0;
		for (Node n : nodesMap.values())
			for (int j = 0; j < n.count; j++)
			{
				real[k] = n.real.intValue();
				predicted[k] = n.predicted.intValue();
				k++;
			}

		int[] classCorrect = new int[numOfClasses], classAll = new int[numOfClasses], classPredictedAll = new int[numOfClasses];

		accuracyTotal.replicaValues = new double[replicas + 1];
		accuracyBalanced.replicaValues = new double[replicas + 1];
		if(mcc == null)mcc = new RangedValue();
		if(auc == null)auc = new RangedValue();
		if(valueByClass == null) valueByClass = new ArrayList<RangedValue>();
		if(sensitivityByClass == null)sensitivityByClass = new ArrayList<RangedValue>();
		
		mcc.replicaValues = new double[replicas + 1];

		for (int j = 0; j < numOfClasses; j++)
		{
			RangedValue value = new RangedValue();
			value.replicaValues = new double[replicas + 1];
			valueByClass.add(value);
		}

		for (int j = 0; j < numOfClasses; j++)
		{
			RangedValue value = new RangedValue();
			value.replicaValues = new double[replicas + 1];
			sensitivityByClass.add(value);
		}

		MersenneTwister random = new MersenneTwister();

		for (int i = 0; i <= replicas; i++)
		{
			for (int j = 0; j < numOfClasses; j++)
				classCorrect[j] = classAll[j] = classPredictedAll[j] = 0;

			for (int j = 0; j < numOfPoints; j++)
			{
				int sample = (i != replicas) ? random.nextInt(numOfPoints) : j; // last result is the real value

				if (real[sample] == predicted[sample])
					classCorrect[real[sample]]++;

				classAll[real[sample]]++;

				if (predicted[sample] < classPredictedAll.length)
					classPredictedAll[predicted[sample]]++;
			}

			accuracyTotal.value = 0D;
			accuracyBalanced.value = 0D;
			for (int j = 0; j < numOfClasses; j++)
			{
				accuracyTotal.value += ((classCorrect[j] * 1.0) / numOfPoints);
				accuracyBalanced.value += ((classCorrect[j] * 1.0) / (classAll[j] * numOfClasses + Double.MIN_VALUE));

				RangedValue v;
				v = sensitivityByClass.get(j);
				v.replicaValues[i] = v.value = ((classCorrect[j] * 1.0) / classAll[j]);

				v = valueByClass.get(j);
				v.replicaValues[i] = v.value = ((classCorrect[j] * 1.0) / classPredictedAll[j]);

			}

			if (numOfClasses == 2)
			{
				tp = classCorrect[0];
				fn = classAll[0] - tp;
				tn = classCorrect[1];
				fp = classAll[1] - tn;

				mcc.value = (1.0 * (tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));
				mcc.value = (mcc.value == 0) ? 0 : (1.0 * (tp * tn - fp * fn) / Math.sqrt(mcc.value)); // to avoid problem with too large
				// 
				mcc.replicaValues[i] = mcc.value;
			}

			accuracyTotal.replicaValues[i] = accuracyTotal.value;
			accuracyBalanced.replicaValues[i] = accuracyBalanced.value;
		}

		for (int j = 0; j < numOfClasses; j++)
		{
			RangedValue v = valueByClass.get(j);
			v.calculateStd();

			v = sensitivityByClass.get(j);
			v.calculateStd();
		}
		accuracyTotal.value *= 100;
		accuracyBalanced.value *= 100;
		accuracyTotal.calculateStd();
		accuracyTotal.std *= 100;
		accuracyBalanced.calculateStd();
		accuracyBalanced.std *= 100;
		if (numOfClasses == 2)
			mcc.calculateStd();
	}

	public void addValue(Double real, Double predicted, Double threshold)
	{
		if(threshold == null) 
			getNode((int)Math.round(real), (int)Math.round(predicted)).count++;
		else // only for two class classification
			getNode((int)Math.round(real), predicted <= threshold?0:1).count++;
	}

	/*
	public void addValue(Integer real, Integer predicted)
	{
		if (predicted < 0)
			predicted = 0;

		if (real < 0)
			real = 0;

		getNode(real, predicted).count++;
	}
	 */
	public long getMaxClass()
	{
		long max = 0;
		for (Node n : nodes)
			if (n.real > max)
				max = n.real;

		return max;
	}

	public static class RangedValue implements Serializable
	{
		private static final long serialVersionUID = 1L;

		Double value = 0D;
		Double std = 0D;

		@XmlTransient
		double[] replicaValues;

		public RangedValue()
		{

		}

		public RangedValue(double value, double std)
		{
			this.value = value;
			this.std = std;
		}

		void calculateValue()
		{
			for (double replica : replicaValues)
				value += replica;
			value /= replicaValues.length;
		}

		void calculateStd()
		{
			std = SetStatistics.prob68Intervals(replicaValues);
		}

		public double getValue()
		{
			return value;
		}

		@XmlAttribute(name = "formatted-std")
		public double getStd()
		{
			if (std == null)
				std = 0.0;
			return Double.parseDouble(NumericalValueStandardizer.getSignificantDigitsStr(std, 1));
		}

		@XmlAttribute(name = "formatted-value")
		public String getFormattedValue()
		{
			if (std == null)
				std = 0.0;

			System.out.println(""+value+" "+std+" --> "+value);
			if(Double.isNaN(value))
				System.out.println(""+value+" "+std+" -*-> "+NumericalValueStandardizer.getWithConfidenceIntervals(value, std));
			return NumericalValueStandardizer.getWithConfidenceIntervals(value, std);
		}
	}

	static class Node implements Serializable
	{
		private static final long serialVersionUID = 4469507048516935373L;

		@XmlAttribute
		public Long real;

		@XmlAttribute
		public Long predicted;

		@XmlAttribute
		public Long count = 0L;

		public String toString()
		{
			return "" + predicted + " for " + real + ": " + count;
		}
	}

	public String toString()
	{
		return "ACC=" + accuracyTotal.getValue() + "+-" + accuracyTotal.getStd() + " " + "BCC=" + accuracyBalanced.getValue() + "+-"
				+ accuracyBalanced.getStd() + " " + "MCC=" + mcc.getValue() + "+-" + mcc.getStd() + " " + " " + "TP=" + tp + " FP=" + fp + " TN=" + tn + " FN="
				+ fn + " " + "N=" + (tp + fp + tn + fn);
	}

	/**
	 *  Prepare a template for the classification report 
	 * @param maxRealValue
	 */

	public void prepareConfusionMatrix(long maxRealValue)
	{

		for (int real = 0; real < maxRealValue; real++)
			for (int predicted = 0; predicted < maxRealValue; predicted++)
				getNode(real, predicted); // all counts will be 0
	}

}