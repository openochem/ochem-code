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

import qspr.entities.ModelMapping;
import qspr.metaserver.configurations.ADConfiguration;
import qspr.modelling.PointStatistics;
import qspr.modelling.SetStatistics;

public class SlidingWindowADFactory extends ADFactory
{
	public SlidingWindowADFactory(SetStatistics statistics, SetStatistics dmSource, ModelMapping mm,
			String dm, ADConfiguration initialConf) {
		super(statistics, dmSource, mm, dm, initialConf);
	}

	public ADConfiguration getAD()
	{
		List<Integer> order = statsWithDM.getDMOrder(dm);
		final int dmIndex = statsWithDM.distancesToModel.indexOf(dm);
		// Average accuracy over a sliding window
		int wndSize = (int)Math.round(order.size() * initialConf.windowSizeInPercent);
		adConfiguration.pointsPerBlock = wndSize;
		Accuracy accuracy = mm.property.isNumeric() ? new RMSE() : new ClassificationAccuracy();
		double sumDm = 0;
		int left = 0;
		for (int i = 0; i < wndSize; i++)
		{
			PointStatistics ps = statsWithPredictions.points.get(order.get(i));
			accuracy.addValue(ps.getRealValue(considerPredicates), ps.predicted);
			sumDm += statsWithDM.points.get(order.get(i)).distancesToModel.get(dmIndex);
		}
		
		// Slide a window...
		while (left + wndSize < order.size())
		{
			int middlePointNum = order.get(Double.valueOf(left + 0.5 * wndSize).intValue());
			adConfiguration.intervals.add(1d * sumDm / wndSize);
			adConfiguration.percents.add(1d * Math.round(1d * (left + wndSize / 2) / order.size() * 10000) / 100);
			adConfiguration.epIds.add(statsWithPredictions.points.get(middlePointNum).id);
			adConfiguration.errors.add(accuracy.getAverageValue());

			PointStatistics ps1 = statsWithPredictions.points.get(order.get(left));
			PointStatistics ps2 = statsWithPredictions.points.get(order.get(left + wndSize));
			PointStatistics dmps1 = statsWithDM.points.get(order.get(left));
			PointStatistics dmps2 = statsWithDM.points.get(order.get(left + wndSize));
			accuracy.removeValue(ps1.getRealValue(considerPredicates), ps1.predicted);
			accuracy.addValue(ps2.getRealValue(considerPredicates), ps2.predicted);
			sumDm -= dmps1.distancesToModel.get(dmIndex);
			sumDm += dmps2.distancesToModel.get(dmIndex);
			
			if (Double.isInfinite(sumDm))
				throw new RuntimeException("Infinite value of DM: " + statsWithPredictions.distancesToModel + " = " + dmps2.distancesToModel);
			
			left++;
		}
		
		return adConfiguration;
	}	
}
