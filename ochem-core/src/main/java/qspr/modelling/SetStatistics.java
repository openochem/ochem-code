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

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.nio.ByteBuffer;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.commons.codec.binary.Hex;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ExperimentalProperty;
import qspr.entities.ModelMapping;
import qspr.interfaces.Descriptable;
import qspr.metaserver.configurations.ADConfiguration;
import qspr.metaserver.protocol.Task;
import qspr.workflow.utils.QSPRConstants;
import qspr.modelling.ClasificationSummary.RangedValue;
import qspr.modelling.ad.ADFactory;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.NumericalValueStandardizer;

@XmlRootElement(name = "set")
public class SetStatistics implements Serializable, Descriptable {
	public static final long serialVersionUID = 1;

	static final int BOOTSTRAP_REPLICAS = 1000;
	static final int BOOTSTRAP_POINT = 10000000;

	@XmlAttribute
	public String setId;

	@XmlTransient
	public double rmse;
	@XmlTransient
	public double r2;
	@XmlTransient
	public double q2;
	@XmlTransient
	public double mae;
	@XmlTransient
	public double average;

	@XmlTransient
	public double rmseStd;
	@XmlTransient
	public double r2Std;
	@XmlTransient
	public double q2Std;
	@XmlTransient
	public double maeStd;
	@XmlTransient
	public double averageStd;

	@XmlTransient
	public int n; // number of samples

	@XmlAttribute
	public int matchedPoints;

	public int errorCount = 0;
	public boolean hasPredicates;

	// For binary classification models
	public ClasificationSummary classificationSummary;
	public Long basketId;

	@XmlElement(name = "dm")
	public List<String> distancesToModel = new ArrayList<String>();

	@XmlElement
	public ADConfiguration adConfiguration;

	@XmlElement(name = "point")
	public List<PointStatistics> points = new ArrayList<PointStatistics>();

	@XmlTransient
	public HashMap<String, ADConfiguration> adConfigurations = new HashMap<String, ADConfiguration>();

	/**
	 * Allows to customize the bootstrapping options for this particular statistics
	 */
	public transient Integer bootstrapReplicas;

	public transient boolean doNotCalculateAUC = false;

	public SetStatistics() {
	}

	public int getNumberOfClasses() {
		return classificationSummary == null?0:classificationSummary.getNumberOfClasses();
	}
	
	public double getPerformanceError() {
		if(errorCount == points.size()) return Double.MAX_VALUE;
		if (classificationSummary == null) return rmse;
		if(classificationSummary.auc == null)classificationSummary.auc = getAUC();
		return 1. - classificationSummary.auc.getValue();
	}

	public boolean hasVirtual() {
		for (PointStatistics ps : points) 
			if(ps.virtual != null) return true; 
		return false;
	}

	public void setPredictions(DataTable dtPredictedValues, ModelMapping modelMapping, Set<Long> recordsInTheSet) throws Exception {
		int index = 0, missed = 0;
		int index_array = -1;
		String col = "";

		distancesToModel.clear();
		// Check weather its multi-property learning
		if (modelMapping.model.modelMappings.size() > 1) 
		{
			col = "" + modelMapping._class;
			missed = Integer.parseInt(col);
			index = dtPredictedValues.getColumnIndex(QSPRConstants.PREDICTION_RESULT_COLUMN + col);
			if(index == -1)throw new UserFriendlyException("\n\nUnexpected failure - can't find predictions for column: " + QSPRConstants.PREDICTION_RESULT_COLUMN + col + "\n\n");
			index_array = dtPredictedValues.getColumnIndex(QSPRConstants.INDIVIDUAL_PREDICTIONS + col);
		} else
			index_array = dtPredictedValues.getColumnIndex(QSPRConstants.INDIVIDUAL_PREDICTIONS);

		String dm = "DM" + col + ":";
		for (String column : dtPredictedValues.getColumns()) {
			if (column.startsWith(dm))
				distancesToModel.add(column.substring(dm.length()));
		}

		if(points.isEmpty())return;

		boolean skipVirtual = hasVirtual();

		List<PointStatistics> virtualpoints = new ArrayList<PointStatistics>();
		HashSet<Long> recordsAdded = new HashSet<Long>();

		if(recordsInTheSet == null || recordsInTheSet.size() == 0) skipVirtual = true;

		for (int record = 0; record <= points.size() && dtPredictedValues.hasMoreRows(); record++) 
		{
			PointStatistics ps = record < points.size()? points.get(record) : null; // last run is required to analyze all data; otherwise we will finish before getting the last records from dtPredictedValues

			if(ps != null && ps.virtual != null) continue;  // virtual points

			Long rowId = -1l;

			do{
				dtPredictedValues.forceNextRow();
				AbstractDataRow row = dtPredictedValues.getCurrentRow();
				rowId = (Long) row.getAttachment(QSPRConstants.RECORD_ID_ATTACHMENT);

				if(rowId == null) {
					if(row.getAttachment(QSPRConstants.CACHED) != null) 
						throw new UserFriendlyException("The data for this basket where uploaded from cache and cannot be used (RECORDID is not available)."+
								" Recalculate basket by disabling cached predictions to add results as the validation set.");
					else
						continue;
				}

				PointStatistics pps = null;

				if(ps == null || !rowId.equals(ps.id) ) {  // Virtual data points - we add virtual predictions for all points that are in the set (but not in others sets)
					if(skipVirtual) continue;
					Double implicit[] = (Double[]) row.getAttachment(QSPRConstants.IMPLICIT_ATTACHMENT);
					if(implicit == null || implicit[missed] == null) continue; // this is missed record
					Long recordId  = (Long) row.getAttachment(QSPRConstants.RECORD_ID_ATTACHMENT);
					if(!recordsInTheSet.contains(recordId) || recordsAdded.contains(recordId))continue; // not from our set or already added
					ExperimentalProperty ep = Repository.record.getRecord(recordId);
					recordsAdded.add(recordId);
					PointStatistics point = new PointStatistics();
					point.predicted = NumericalValueStandardizer.getSignificantDigitsDouble(((Double) row.getValue(index)).doubleValue(),NumericalValueStandardizer.SIGNIFICANT_DIGITS);
					point.real = implicit[missed];
					point.virtual = true;
					point.molWeight = ep.molecule != null?  ep.molecule.getMolWeight() : 0.;
					ExperimentalProperty pp = new ExperimentalProperty(ep.molecule);
					pp.property = Repository.property.getPropertyById(modelMapping.property.id);
					point.id = -ep.id;
					point.moleculeId = ep.molecule.mapping1.id;
					virtualpoints.add(point);
					pps = point; // working with a virtual record
				}
				else {
					pps = ps; // found our record!
					if(!skipVirtual) recordsAdded.add(pps.id); // to avoid duplicates if this record in the validation set 
				}

				if (Task.ERROR.equals(row.status))
					pps.error = row.detailedStatus;
				else 
				{
					if (row.getValue(index) == null)
						throw new UserFriendlyException("Null prediction value in row "	+ dtPredictedValues.currentRow + "!");
					pps.predicted = NumericalValueStandardizer.getSignificantDigitsDouble(((Double) row.getValue(index)).doubleValue(),NumericalValueStandardizer.SIGNIFICANT_DIGITS);

					if (index_array != -1)
						pps.ensemblePredictions = (float[])row.getValue(index_array);

					pps.error = null;
					Iterator<String> iterDms = distancesToModel.iterator();
					while (iterDms.hasNext()) 
					{
						String curDM = iterDms.next();
						if (dtPredictedValues.getValue(dm + curDM) != null)
							pps.distancesToModel.add(((Double) dtPredictedValues.getValue(dm + curDM)).doubleValue());
						else
							iterDms.remove();
					}
				}
			}
			while( (ps == null || rowId == null || !rowId.equals(ps.id)) && dtPredictedValues.hasMoreRows()); 

			if(ps != null && rowId != null && !rowId.equals(ps.id)) 
				throw new UserFriendlyException("Some data were lost: please, report this bug");
		}

		if(!skipVirtual){ // adding missed 

		}

		points.addAll(virtualpoints);

		recalculateStatistics(modelMapping);
	}

	public SetStatistics(Basket basket, ModelMapping modelMapping) {
		basketId = basket.id;
		for (BasketEntry entry : basket.entries)
			if (modelMapping.matches(entry.ep))
				points.add(new PointStatistics(entry.ep, modelMapping));
	}

	public SetStatistics(SetStatistics set1,
			SetStatistics set2) {
		if(set1!=null) {points.addAll(set1.points); basketId = set1.basketId;}
		if(set2!=null) {points.addAll(set2.points); basketId = set2.basketId;}
	}

	public void recalculateStatistics(ModelMapping modelMapping) {
		recalculateStatistics(modelMapping, null);
	}

	private void statisticsClassification(ModelMapping modelMapping, PointSelector selector)
	{
		classificationSummary = new ClasificationSummary();
		long maxRealValue = 0;

		if (modelMapping.property.isQualitative())
			for (PointStatistics ps : points)
				if (ps.real > maxRealValue)
					maxRealValue = Math.round(ps.real);

		hasPredicates = false;
		n = 0;
		matchedPoints = 0;

		//When can maxRealValue be 0?
		maxRealValue = (maxRealValue==0) ? 2 : maxRealValue+1;

		if(!Globals.considerPredicates() || maxRealValue > 2) { // for multiclass no change!
			modelMapping.classificationThreshlod = null;
		}
		else
			if(modelMapping.classificationThreshlod == null) 
				modelMapping.classificationThreshlod = selectThreshold(modelMapping,selector);

		classificationSummary.prepareConfusionMatrix(maxRealValue);
		errorCount = 0;
		for (PointStatistics ps : points) 
		{
			n++;

			if (ps.error != null)
			{
				errorCount++;
				continue;
			}

			if (selector != null && !selector.pointMatches(ps, this))
				continue;

			classificationSummary.addValue(ps.real < 0?0:ps.real, ps.predicted < 0 && modelMapping.classificationThreshlod == null ? 0:ps.predicted, 
					modelMapping.classificationThreshlod);

			matchedPoints++;
		}

		classificationSummary.calculateStatistics(getReplicasCount());

		if (!doNotCalculateAUC)
			classificationSummary.auc = getAUC();
		//
	}

	private Double selectThreshold(ModelMapping modelMapping, PointSelector selector) {
		Map <Double,Integer>negatives = new HashMap<Double,Integer>();
		Map <Double,Integer>positives = new HashMap<Double,Integer>();
		Set <Double> thresholds = new TreeSet<Double>();
		int negative=0,positive=0;

		for (PointStatistics ps : points)
			if (ps.error == null)
			{
				if (selector != null && !selector.pointMatches(ps, this))
					continue;
				double predicted = ps.predicted;

				if(!thresholds.contains(predicted)) {
					negatives.put(predicted, 0);
					positives.put(predicted, 0);
					thresholds.add(predicted);
				}

				if(ps.real < 0.5){
					negative++;
					negatives.put(predicted, negatives.get(predicted)+1);
				}else {
					positive++;
					positives.put(predicted, positives.get(predicted)+1);
				}
			}

		double bestThreshold = 0, bestAccuracy = 0;
		double startPositives = positive, startNegatives = 0;

		for(Double threshold:thresholds) {
			startPositives -= positives.get(threshold);
			startNegatives += negatives.get(threshold);
			double balanced = (startPositives/positive + startNegatives/negative)/2;
			if(balanced > bestAccuracy) {
				bestAccuracy = balanced;
				bestThreshold = threshold;
			}
		}

		System.out.println("best threshold " + bestThreshold + " bestAccuracy = " + bestAccuracy);

		return bestThreshold;
	}

	public RangedValue getAUC()
	{
		long timer = System.nanoTime();
		RangedValue initial = getAUC(1);
		if(getReplicasCount() <= 1) return initial;
		RangedValue rv = getAUC(getReplicasCount());
		rv.value = initial.value;
		logger.info("AUC resampled from "+ getReplicasCount() +" curves on "+points.size()+" points in "+(System.nanoTime() - timer)/1E6+" ms " + rv.getFormattedValue());
		return rv;
	}


	private RangedValue getAUC(int numResamples) {
		RangedValue rv = new RangedValue();
		rv.replicaValues = new double[numResamples];
		for (int i=0; i<numResamples; i++)
		{
			ROCCurve curve = getROCCurve(i); // first curve with seed 0 - no resampling
			ROCPoint p = null;
			double area = 0D; 
			for (ROCPoint point : curve.points) 
			{
				if (p != null)
					area += (point.x - p.x)*(p.y + (point.y - p.y) / 2);
				p = point;
			}
			rv.replicaValues[i] = area;
		}

		rv.calculateValue();
		rv.calculateStd();
		return rv;
	}

	public ROCCurve getROCCurve()
	{
		return getROCCurve(-1);
	}

	private ROCCurve getROCCurve(int resampleSeed)
	{
		ROCCurve curve = new ROCCurve(setId);

		List<PointStatistics> l = new ArrayList<PointStatistics>();

		Random r = null;
		if (resampleSeed > 0)
			r = new Random(resampleSeed);

		for (int i=0; i < points.size(); i++)
		{
			PointStatistics ps = null;
			if (r != null)
				ps = points.get(r.nextInt(points.size()));
			else
				ps = points.get(i);

			if (ps.error == null)
				l.add(ps);
		}

		Comparator<PointStatistics> c = new Comparator<PointStatistics>()
		{
			@Override
			public int compare(PointStatistics a1, PointStatistics a2) 
			{
				return Double.valueOf(a1.predicted).compareTo(a2.predicted);
			}
		};
		Collections.sort(l, c);
		int positive = 0, negative = 0;
		int tp = 0, fn = 0;

		for (PointStatistics ps : l)
			if (Double.valueOf(ps.real).intValue() == 1)
				positive++;
			else
				negative++;

		double predicted = Double.NaN;
		curve.addPoint(0D, 0D, -1L);
		for (PointStatistics ps : l) 
		{
			if (Double.valueOf(ps.real).intValue() == 1) //positive
				tp++;
			else
				fn++;

			if (Double.isNaN(predicted) || Math.abs(predicted - ps.predicted) > 1E-6)
				curve.addPoint(tp * 1.0D / positive, fn * 1.0D / negative, ps.id);
			else
				curve.updatePoint(tp * 1.0D / positive, fn * 1.0D / negative, ps.id);

			predicted = ps.predicted;
		}
		return curve;
	}

	private void statisticsRegression(ModelMapping modelMapping, PointSelector selector)
	{
		classificationSummary = null;
		int noOfSamples = 0;
		boolean considerPredicates = Globals.considerPredicates();
		double real[] = new double[points.size()], predicted[] = new double[points.size()];
		hasPredicates = false;
		n = 0;
		matchedPoints = 0;
		errorCount = 0;
		Set<Long> molecules = new HashSet<Long>();
		for (PointStatistics ps : points) 
		{
			n++;

			if (ps.error != null)
			{
				errorCount++;
				continue;
			}

			if (selector != null && !selector.pointMatches(ps, this))
				continue;

			real[noOfSamples] = ps.getRealValue(considerPredicates);
			predicted[noOfSamples++] = ps.predicted;

			if (ps.predicate != null && !ps.predicate.equals("="))
				hasPredicates = true;

			matchedPoints++;

			molecules.add(ps.moleculeId);
		}

		long start = Calendar.getInstance().getTimeInMillis();
		bootstrap(real, predicted, noOfSamples, getReplicasCount());
		logger.info("Statistics bootstrap with "+getReplicasCount() + " calculated in " + (Calendar.getInstance().getTimeInMillis() - start)+"ms");
	}

	private int getReplicasCount()
	{
		int replicas = bootstrapReplicas != null?bootstrapReplicas:BOOTSTRAP_REPLICAS; // maximum number of replicas
		int maxreplicas =  points.size() > 0 ?BOOTSTRAP_POINT / points.size():2;
		if(maxreplicas < 2) maxreplicas = 2;
		return replicas < maxreplicas ? replicas : maxreplicas;
	}


	public void recalculateStatistics(ModelMapping modelMapping, PointSelector selector) 
	{
		if (modelMapping.property.isQualitative())
			statisticsClassification(modelMapping, selector);
		else
			statisticsRegression(modelMapping, selector);
	}

	/**
	 * Estimates variance of predictions using bootstrap
	 */

	private void bootstrap(double[] real, double[] predicted, int noOfSamples,
			int replicas) {

		if(noOfSamples<1)return;

		double r2Boot[] = new double[replicas + 1], 
				rmseBoot[] = new double[replicas + 1], 
				q2Boot[] = new double[replicas + 1], 
				maeBoot[] = new double[replicas + 1],
				averageBoot[] = new double[replicas + 1];

		for (int i = 0; i <= replicas; i++) {
			rmse = mae = average = 0;
			double x_sum = 0, y_sum = 0, xx_sum = 0, yy_sum = 0, xy_sum = 0, error, realVal, predictedVal;

			for (int j = 0; j < noOfSamples; j++) {
				int sample = i != replicas ? (int) (Math.random() * noOfSamples)
						: j; // last result is the real value
				realVal = real[sample];
				predictedVal = predicted[sample];
				error = realVal - predictedVal;
				rmse += error * error;
				mae += error > 0 ? error : -error;
				average += error;
				x_sum += realVal;
				y_sum += predictedVal;
				xx_sum += realVal * realVal;
				yy_sum += predictedVal * predictedVal;
				xy_sum += realVal * predictedVal;
			}

			rmse = Math.sqrt(rmse / noOfSamples);
			mae /= noOfSamples;
			average /= noOfSamples;
			r2 = ((noOfSamples * xy_sum - x_sum * y_sum) / (Math
					.sqrt((noOfSamples * xx_sum - x_sum * x_sum +0.)
							* (noOfSamples * yy_sum - y_sum * y_sum +0.))));
			q2 = 1 - (yy_sum + xx_sum - 2 * xy_sum)
					/ (xx_sum - x_sum / noOfSamples * x_sum);
			r2 = r2 * r2;
			r2 = r2 < -1 ? -1 : r2 > 1 ? 1 : r2; // rounding to avoid erroneous
			// values due to rounding
			// problems
			q2 = q2 < 0 ? 0 : q2 > 1 ? 1 : q2;

			r2Boot[i] = r2;
			q2Boot[i] = q2;
			rmseBoot[i] = rmse;
			maeBoot[i] = mae;
			averageBoot[i] = average;
		}

		r2Std = prob68Intervals(r2Boot);
		q2Std = prob68Intervals(q2Boot);
		rmseStd = prob68Intervals(rmseBoot);
		maeStd = prob68Intervals(maeBoot);
		averageStd = prob68Intervals(averageBoot);
	}

	static double prob68Intervals(double[] bootstrap) {
		if(bootstrap.length <= 1)
			return 0; // just some value for cases when we do not do boostrapping

		Arrays.sort(bootstrap);
		int n1 = (int) Math.round(0.84 * bootstrap.length); n1 = n1 < bootstrap.length ? n1 : bootstrap.length - 1;
		int n2 = (int) Math.round(0.16 * bootstrap.length);
		return (bootstrap[n1] - bootstrap[n2]) / 2.;
	}

	// Get a sample with replacement (usually for a significancy test)
	public SetStatistics getSample(List<Integer> sample) {
		SetStatistics sampleStat = new SetStatistics();
		for (int i = 0; i < sample.size(); i++) {
			sampleStat.points.add(points.get(sample.get(i)));
		}

		return sampleStat;
	}

	public SetStatistics setId(String id) {
		setId = id;
		return this;
	}

	public double getDM(int pointNum, String dm) {
		return points.get(pointNum).distancesToModel.get(dm == null ? 0
				: distancesToModel.indexOf(dm));
	}

	public void setDM(int pointNum, String dm, double value) {
		int dmIndex = distancesToModel.indexOf(dm);
		if (dmIndex == -1) {
			distancesToModel.add(dm);
			dmIndex = distancesToModel.size() - 1;
		}

		PointStatistics point = points.get(pointNum);
		while (point.distancesToModel.size() - 1 < dmIndex)
			point.distancesToModel.add(new Double(0));
		point.distancesToModel.set(dmIndex, value);
	}

	public Map<String, Object> getParameters() {
		Map<String, Object> parameters = new HashMap<String, Object>();
		parameters.put("RMSE", rmse);
		parameters.put("r2", r2);
		return parameters;
	}

	public ADConfiguration getADConfiguration(String dm,
			SetStatistics dmSource, ModelMapping mm,
			ADConfiguration configuration) {
		return ADFactory.getInstance(this, dmSource, mm, dm, configuration)
				.getAD();
	}

	private Map<String, List<Integer>> orders;

	public List<Integer> getDMOrder(String dm) {
		if (orders == null)
			orders = new HashMap<String, List<Integer>>();
		List<Integer> order = orders.get(dm);
		if (order != null)
			return order;

		order = new ArrayList<Integer>();
		final int dmIndex = distancesToModel.indexOf(dm);
		int i = 0;
		for (PointStatistics ps : points) {
			if (ps.error == null)
				order.add(i);
			i++;
		}

		// Calculate the DM-based order
		Collections.sort(order, new Comparator<Integer>() {
			public int compare(Integer o1, Integer o2) {
				Double dm1 = points.get(o1).distancesToModel.get(dmIndex);
				Double dm2 = points.get(o2).distancesToModel.get(dmIndex);
				if (dm1.equals(dm2))
					return 0;
				else
					return dm1 > dm2 ? 1 : -1;
			}
		});
		orders.put(dm, order);
		return order;
	}

	// Reorder the points according to points of another SetStatistics
	public void reorderAccordingTo(SetStatistics template) throws Exception {
		List<PointStatistics> tmpPoints = new ArrayList<PointStatistics>();
		Map<Long, Integer> rowById = new HashMap<Long, Integer>();
		for (int i = 0; i < points.size(); i++)
			rowById.put(points.get(i).id, i);

		for (int i = 0; i < template.points.size(); i++) {
			Integer row = rowById.get(template.points.get(i).id);
			if (row == null)
				throw new Exception("Selected DM is incompatible with the data");
			tmpPoints.add(points.get(row));
		}

		points = tmpPoints;
		orders = null;
	}

	public String toString() {
		if (classificationSummary != null)
			return classificationSummary.toString();
		else
			return "RMSE: " + getRmseStr() + ", R2: " + getR2Str() + ", MAE:"
			+ getMaeStr();
	}

	@XmlElement(name = "r2")
	public String getR2Str() {
		if (r2 < 0.01)
			return "0";
		return NumericalValueStandardizer.getWithConfidenceIntervals(r2, r2Std);
	}

	protected void setR2Str(String val) {
		r2 = Double.valueOf(val);
	}

	@XmlElement(name = "q2")
	public String getQ2Str() {
		if (q2 < 0.01)
			return "0";
		return NumericalValueStandardizer.getWithConfidenceIntervals(q2, q2Std);
	}

	protected void setQ2Str(String val) {
		q2 = Double.valueOf(val);
	}

	@XmlElement(name = "rmse")
	public String getRmseStr() {
		return NumericalValueStandardizer.getWithConfidenceIntervals(rmse,
				rmseStd);
	}

	protected void setRmseStr(String val) {
		rmse = Double.valueOf(val);
	}

	@XmlElement(name = "mae")
	public String getMaeStr() {
		return NumericalValueStandardizer.getWithConfidenceIntervals(mae,
				maeStd);
	}

	protected void setMaeStr(String val) {
		mae = Double.valueOf(val);
	}

	@XmlElement(name = "average")
	public String getAverageStr() {
		return NumericalValueStandardizer.getWithConfidenceIntervals(average,
				averageStd);
	}

	protected void setAverageStr(String val) {
		average = Double.valueOf(val);
	}

	@XmlElement(name = "r2std")
	protected String getR2StdStr() {
		return NumericalValueStandardizer.getSignificantDigitsStr(r2Std, 1);
	}

	protected void setR2StdStr(String val) {
		r2Std = Double.valueOf(val);
	}

	@XmlElement(name = "q2std")
	protected String getQ2StdStr() {
		return NumericalValueStandardizer.getSignificantDigitsStr(q2Std, 1);
	}

	public void setQ2StdStr(String val) {
		q2Std = Double.valueOf(val);
	}

	@XmlElement(name = "rmsestd")
	protected String getRmseStdStr() {
		return NumericalValueStandardizer.getSignificantDigitsStr(rmseStd, 1);
	}

	protected void setRmseStdStr(String val) {
		rmseStd = Double.valueOf(val);
	}

	@XmlElement(name = "maestd")
	protected String getMaeStdStr() {
		return NumericalValueStandardizer.getSignificantDigitsStr(maeStd, 1);
	}

	protected void setMaeStdStr(String val) {
		maeStd = Double.valueOf(val);
	}

	@XmlAttribute(name = "size")
	protected int getRowsSize() {
		if (points == null)
			return 0;
		return points.size() - errorCount;
	}

	@XmlElement(name = "accuracy")
	public String getAccuracyStr() {
		if (classificationSummary != null)
			return classificationSummary.accuracyTotal.getFormattedValue();
		return null;
	}

	@XmlElement(name = "balancedAccuracy")
	public String getBalancedAccuracyStr() {
		if (classificationSummary != null)
			return classificationSummary.accuracyBalanced.getFormattedValue();
		return null;
	}

	@XmlElement(name = "auc")
	public String getAucAccuracyStr() {
		if (classificationSummary != null)
			return classificationSummary.auc.getFormattedValue();
		return null;
	}

	@XmlElement(name = "mcc")
	public String getMCCStr() {
		if (classificationSummary != null)
			return classificationSummary.mcc.getFormattedValue();
		return null;
	}


	/**
	 * Get the MD5 hash of the set based on the point IDs.
	 * Allows to identify changes in the sets.
	 */
	@XmlTransient
	public String getHash() throws NoSuchAlgorithmException{
		MessageDigest md5 = MessageDigest.getInstance("MD5");
		md5.reset();
		for (PointStatistics ps : points)
			if(ps.virtual == null)
				md5.update(ByteBuffer.allocate(8).putLong(ps.id).array());
		return String.valueOf(Hex.encodeHex(md5.digest()));
	}

	public void saveToFile(String fileName) throws FileNotFoundException {
		PrintWriter pw = new PrintWriter(new File(fileName));
		pw.println("SMILES,VALUE");
		for (PointStatistics ps : points) {
			ExperimentalProperty ep = (ExperimentalProperty) Globals.session()
					.get(ExperimentalProperty.class, ps.id);
			pw.println(String.format("\"%s\",\"%f\"", ep.molecule.getData()
					.replaceAll("\n", "newline"), ps.real));
		}

		pw.flush();
		pw.close();
	}

	@XmlTransient
	public List<Long> getEpIds() {
		List<Long> res = new ArrayList<Long>();
		for (int i = 0; i < points.size(); i++)
			res.add(points.get(i).id);
		return res;
	}

	// TODO: Consider moving statistics to database tables - SetStatistics and
	// PointStatistics
	// Thus - more transparent structure (avoid serializables) and possibility
	// to make conditional queries to statistics
	// Midnighter

	private static transient Logger logger = LogManager.getLogger(SetStatistics.class);

}
