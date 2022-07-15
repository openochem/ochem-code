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

package com.eadmet.business;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;

import qspr.modelling.PointStatistics;

public class PredictionNeighbor implements Comparable<PredictionNeighbor>
{
	public PointStatistics trainingSetPoint;
	public Double distance;
	
	private PredictionNeighbor()
	{
		
	}
	
	public PredictionNeighbor(PointStatistics ps, double distance)
	{
		this.trainingSetPoint = ps;
		this.distance = distance;
	}
	
	public static PredictionNeighbor getNeighbour(float[] predictionVector, PointStatistics trainingSetPoint, String type)
	{
		PredictionNeighbor de = new PredictionNeighbor();
		de.trainingSetPoint = trainingSetPoint;
		
		DistanceCalculator dc;
		if (type.equals("correlation"))
			dc = new CorrelDistanceCalculator();
		else if (type.equals("rank-correlation"))
			dc = new RankCorrelDistanceCalculator();
		else //"euclidean"
			dc = new EuclideanDistanceCalculator();
		
		if (predictionVector == null || trainingSetPoint.ensemblePredictions == null)
			return null;
		
		int validLength = 0;
		
		for (int i=0; i<Math.min(predictionVector.length, trainingSetPoint.ensemblePredictions.length); i++)
			if (!Float.isNaN(predictionVector[i]) && !Float.isNaN(trainingSetPoint.ensemblePredictions[i]))
				validLength++;
		
		if (validLength < 3)
			return null;
		
		double[] src = new double[validLength];
		double[] dst = new double[validLength];
		
		for (int i=0; i<Math.min(predictionVector.length, trainingSetPoint.ensemblePredictions.length); i++)
			if (!Float.isNaN(predictionVector[i]) && !Float.isNaN(trainingSetPoint.ensemblePredictions[i]))
			{
				validLength--;
				src[validLength] = predictionVector[i];
				dst[validLength] = trainingSetPoint.ensemblePredictions[i];
			}

		de.distance = dc.getDistance(src, dst);
		return de;
	}

	@Override
	public int compareTo(PredictionNeighbor arg) 
	{
		return -distance.compareTo(arg.distance);
	}
}

abstract class DistanceCalculator
{
	public abstract double getDistance(double[] vec1, double[] vec2);
}

class RankCorrelDistanceCalculator extends DistanceCalculator
{
	SpearmansCorrelation sc = new SpearmansCorrelation();
	
	@Override
	public double getDistance(double[] vec1, double[] vec2)
	{
		return sc.correlation(vec1, vec2);
	}
}

class CorrelDistanceCalculator extends DistanceCalculator
{
	PearsonsCorrelation pc = new PearsonsCorrelation();
	
	@Override
	public double getDistance(double[] vec1, double[] vec2)
	{
		return pc.correlation(vec1, vec2);
	}
}

class EuclideanDistanceCalculator extends DistanceCalculator
{
	@Override
	public double getDistance(double[] vec1, double[] vec2)
	{
		double sum = 0;
		for (int i=0; i<vec1.length; i++)
			sum += (vec1[i] - vec2[i])*(vec1[i] - vec2[i]);
		return -Math.sqrt(sum) / vec1.length; //Bigger = better, so that compatible to Correl
	}
}

