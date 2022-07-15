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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.workflow.utils.QSPRConstants;

/**
 * A selector for the best predictions according to a DM
 * @author midnighter
 *
 */
public class BestPredictionsPointSelector implements PointSelector
{
	private Map<SetStatistics, Double> dmThreasolds = new HashMap<SetStatistics, Double>();
	private static final Logger logger = LogManager.getLogger(BestPredictionsPointSelector.class);
	private double percentageOfPredictions;
	private boolean useTrainingSetThreshold;
	private int dmNum;
	private SetStatistics ssTraining;
	private Double fixedDMThreshold;
	
	public boolean pointMatches(PointStatistics ps, SetStatistics ss)
	{
		if (ss.setId.equals(QSPRConstants.TRAINING))
			ssTraining = ss;
		return ps.distancesToModel == null || ps.distancesToModel.isEmpty() || ps.distancesToModel.get(0) <= getDMThreshold(ss);
	}
	
	private double getDMThreshold(SetStatistics ss)
	{
		if (fixedDMThreshold != null)
			return fixedDMThreshold;
		
		if (useTrainingSetThreshold && !ss.equals(ssTraining))
			return getDMThreshold(ssTraining);
		
		Double res = dmThreasolds.get(ss);
		if (res != null)
			return res;
		List<Double> dms = new ArrayList<Double>();
	
		for (PointStatistics ps : ss.points)
			if (ps.distancesToModel != null && !ps.distancesToModel.isEmpty())
				dms.add(ps.distancesToModel.get(dmNum)); // We always consider the first DM. Can be changed later.
		
		Collections.sort(dms);
		
		res = dms.get((int)(1 * Math.round(dms.size() * percentageOfPredictions) - 1));
		dmThreasolds.put(ss, res);
		
		logger.info("DM threshold for best " + Math.round(percentageOfPredictions * 100) + "% in set " + ss.setId + " is " + res);
		
		return res;
	}
	
	public BestPredictionsPointSelector(double percentageOfPredictions, boolean useTrainingSetThreshold, int dmNum)
	{
		this.percentageOfPredictions = percentageOfPredictions;
		this.useTrainingSetThreshold = useTrainingSetThreshold;
		this.dmNum = dmNum;
	}
	
	public BestPredictionsPointSelector(double fixedDMThreshold)
	{
		this.fixedDMThreshold = fixedDMThreshold;
	}
}
