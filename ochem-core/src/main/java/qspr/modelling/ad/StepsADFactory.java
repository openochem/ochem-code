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

package qspr.modelling.ad;

import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.entities.ModelMapping;
import qspr.metaserver.configurations.ADConfiguration;
import qspr.modelling.PointStatistics;
import qspr.modelling.SetStatistics;

// AD, based on "MGD-steps" / Midnighter
public class StepsADFactory extends ADFactory
{
	private static transient final Logger logger = LogManager.getLogger(StepsADFactory.class);
	
	public StepsADFactory(SetStatistics statistics, SetStatistics dmSource, ModelMapping mm,
			String dm, ADConfiguration initialConf) {
		super(statistics, dmSource, mm, dm, initialConf);
	}

	@Override
	public ADConfiguration getAD() 
	{
		List<Integer> order = statsWithDM.getDMOrder(dm);
		int dmIndex = statsWithDM.distancesToModel.indexOf(dm);
		
		Accuracy accuracy = Accuracy.forProperty(mm.property);
		Accuracy prevAccuracy = null;
		for (int i = 0; i < order.size(); i++) 
		{
			PointStatistics ps = statsWithPredictions.points.get(order.get(i));
			PointStatistics psDm = statsWithDM.points.get(order.get(i));
			Double dm = psDm.distancesToModel.get(dmIndex);
			accuracy.addValue(ps.getRealValue(considerPredicates), ps.predicted);
			
			boolean accuracyDecreased = prevAccuracy == null || !accuracy.isBetterThan(prevAccuracy);
			int pointsLeft = order.size() - i - 1;
			if (pointsLeft == 0 || ((accuracy.count >= initialConf.pointsPerBlock && pointsLeft > initialConf.pointsPerBlock / 2) && accuracyDecreased))
			{
				double accValue = accuracyDecreased ? accuracy.getAverageValue() : prevAccuracy.getAverageValue();
				adConfiguration.errors.add(accValue);
				adConfiguration.intervals.add(dm);
				adConfiguration.percents.add(1d * Math.round(1d * i / order.size() * 10000) / 100);
				prevAccuracy = accuracy;
				logger.info("AD step: " + accValue + " upto " + dm + " of " + dm + " ("+accuracy.count+" points)");
				accuracy = Accuracy.forProperty(mm.property);
			}
		}
		return adConfiguration;
	}
}