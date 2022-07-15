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

package qspr.modelling.applier;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.Random;

import javax.xml.bind.annotation.XmlRootElement;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;

import qspr.entities.ModelMapping;
import qspr.workflow.utils.QSPRConstants;
import qspr.metaserver.util.aggregator.AggregatedValue;
import qspr.modelling.ApplicabilityDomain;
import qspr.modelling.ClasificationSummary;
import qspr.modelling.ModelStatistics;
import qspr.modelling.PointStatistics;
import qspr.modelling.SetStatistics;
import qspr.workflow.datatypes.DataTable;

/**
 * The prediction results overview.
 * This class is under construction. Lots of related code is scattered in ModelController
 * 
 * @author midnighter
 *
 */
@XmlRootElement
public class PredictionResults
{
	//public transient DataTable dtPredictions;
	public transient List<PropertyPrediction> predictions = new ArrayList<PropertyPrediction>();
	public int propNum;

	public transient ApplicabilityDomain ad;

	public ModelMapping modelMapping;

	public boolean removeErrors = true;

	/**
	 * The statistics based on error estimates simulated via bootstrapping
	 */
	public SetStatistics simulatedStatistics;

	public synchronized void setPredictions(DataTable dtPredictions, ApplicabilityDomain ad, ModelMapping mm, int propNum)
	{
		int col = dtPredictions.indexOfColumn(QSPRConstants.PREDICTION_RESULT_COLUMN + propNum);
		if (col == -1)
			col = 0;

		dtPredictions.reset();
		while (dtPredictions.nextRow())
		{
			if (dtPredictions.getCurrentRow().isError())
				if (removeErrors)
					continue;
			PropertyPrediction pp = new PropertyPrediction();
			if (dtPredictions.getCurrentRow().isError())
				pp.setError(dtPredictions.getCurrentRow().detailedStatus);
			else
			{
				if (ad != null)
				{
					pp.setAccuracy(ad.getPredictedError(dtPredictions.currentRow));
					pp.setInsideAD(ad.getPredictedDM(dtPredictions.currentRow) < ad.dmThreshold);
				}
				pp.setValue((Double) dtPredictions.getValue(col));
			}
			predictions.add(pp);

		}

		this.propNum = propNum;
		this.modelMapping = mm;
		this.ad = ad;
	}

	public void setPredictions(SetStatistics ss, ApplicabilityDomain ad, ModelMapping mm, int propNum)
	{
		int dmIndex = ad != null ? ss.distancesToModel.indexOf(ad.dmName) : 0;
		for (PointStatistics ps : ss.points)
		{
			if (ps.error != null)
				if (removeErrors)
					continue;
			PropertyPrediction pp = new PropertyPrediction();
			if (ad != null)
				pp.setAccuracy(ad.getPredictedError(ps.distancesToModel.get(dmIndex)));
			pp.setValue(ps.predicted);
			predictions.add(pp);
		}

		this.propNum = propNum;
		this.modelMapping = mm;
		this.ad = ad;
	}

	/**
	 * Simulate the predictive statistics based on AD-based simulated real values
	 */
	
	private static transient Logger logger = LogManager.getLogger(PredictionResults.class);
	
	public void bootstrapStatistics()
	{
		List<SetStatistics> sss = new ArrayList<SetStatistics>();

		long maxClass = 0;
		if (modelMapping.property.isQualitative())
		{
			SetStatistics ss = ModelStatistics.get(modelMapping).sets.get(0);
			if (ss.classificationSummary == null)
				ss.recalculateStatistics(modelMapping);
			maxClass = ss.classificationSummary.getMaxClass();
		}

		int replicas = 1000000 / (predictions.size() == 0?1000000:predictions.size());
		if(replicas > 1000) replicas = 1000; 
		if(replicas < 1) replicas = 1;

		long start = Calendar.getInstance().getTimeInMillis();
		
		for (int i = 0; i < replicas; i++)
		{
			SetStatistics ss = new SetStatistics();
			ss.bootstrapReplicas = QSPRConstants.NO_REPLICAS;
			ss.doNotCalculateAUC = true;

			Random random = new Random();

			//dtPredictions.reset();
			//while (dtPredictions.nextRow())
			for (int k = 0; k < predictions.size(); k++)	
			{
				PropertyPrediction prediction = predictions.get(k);
				if (prediction.getError() != null)
					continue;

				PointStatistics ps = new PointStatistics(prediction.getValue());

				if (modelMapping.property.isQualitative())
				{
					// Makes sense only for classification, amiright?
					if (ps.predicted < 0)
						ps.predicted = 0;
					if (ps.predicted >= maxClass)
						ps.predicted = maxClass;

					// Simulate classification error
					if (random.nextDouble() > prediction.getAccuracy())
						ps.real = maxClass - ps.predicted;
					else
						ps.real = ps.predicted;
				}
				else
					ps.real = random.nextGaussian() * prediction.getAccuracy() + ps.predicted;

				ss.points.add(ps);
			}

			Logger logger = LogManager.getLogger(SetStatistics.class);
			Configurator.setLevel(logger.getName(), Level.ERROR);
			// LogManager.getLogger(SetStatistics.class).removeAllAppenders(); N.B.! 
			ss.recalculateStatistics(modelMapping);
			ss.points = null; // Release memory
			sss.add(ss);
		}

		logger.info("Profiling: Statistics bootstrap with "+replicas + " calculated in " + (Calendar.getInstance().getTimeInMillis() - start)+"ms");

		simulatedStatistics = new SetStatistics();

		AggregatedValue rmse = new AggregatedValue(true);
		AggregatedValue mae = new AggregatedValue(true);
		AggregatedValue r2 = new AggregatedValue(true);
		AggregatedValue q2 = new AggregatedValue(true);
		AggregatedValue accuracy = new AggregatedValue(true);
		AggregatedValue balancedAccuracy = new AggregatedValue(true);
		AggregatedValue mcc = new AggregatedValue(true);

		for (SetStatistics ss : sss)
		{
			simulatedStatistics.matchedPoints = ss.matchedPoints;
			if (!modelMapping.property.isQualitative())
			{
				rmse.addValue(ss.rmse);
				mae.addValue(ss.mae);
				r2.addValue(ss.r2);
				q2.addValue(ss.q2);
			}
			else
			{
				accuracy.addValue(ss.classificationSummary.accuracyTotal.getValue());
				balancedAccuracy.addValue(ss.classificationSummary.accuracyBalanced.getValue());
				mcc.addValue(ss.classificationSummary.mcc.getValue());
			}
		}

		if (!modelMapping.property.isQualitative())
		{
			simulatedStatistics.rmse = rmse.getRMSE();
			simulatedStatistics.rmseStd = rmse.get95ConfidenceInterval();

			simulatedStatistics.mae = mae.getMean();
			simulatedStatistics.maeStd = mae.get95ConfidenceInterval();

			simulatedStatistics.r2 = r2.getMean();
			simulatedStatistics.r2Std = r2.get95ConfidenceInterval();

			simulatedStatistics.q2 = q2.getMean();
			simulatedStatistics.q2Std = q2.get95ConfidenceInterval();
		}
		else
		{
			simulatedStatistics.classificationSummary = new ClasificationSummary();
			simulatedStatistics.classificationSummary.accuracyTotal = new ClasificationSummary.RangedValue(accuracy.getMean(), accuracy.get95ConfidenceInterval());
			simulatedStatistics.classificationSummary.accuracyBalanced = new ClasificationSummary.RangedValue(balancedAccuracy.getMean(), balancedAccuracy.get95ConfidenceInterval());
			simulatedStatistics.classificationSummary.mcc = new ClasificationSummary.RangedValue(mcc.getMean(), mcc.get95ConfidenceInterval());
		}
	}
}
