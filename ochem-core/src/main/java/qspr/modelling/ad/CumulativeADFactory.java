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

public class CumulativeADFactory extends ADFactory
{
	public CumulativeADFactory(SetStatistics statistics, SetStatistics dmSource, ModelMapping mm,
			String dm, ADConfiguration initialConf) {
		super(statistics, dmSource, mm, dm, initialConf);
	}

	public ADConfiguration getAD()
	{
		List<Integer> order = statsWithDM.getDMOrder(dm);
		final int dmIndex = statsWithDM.distancesToModel.indexOf(dm);
		Double prevDm = statsWithDM.points.get(order.get(0)).distancesToModel.get(dmIndex);
		Accuracy accuracy = Accuracy.forProperty(mm.property);
		for (int i = 0; i < order.size(); i++) 
		{
			PointStatistics ps = statsWithPredictions.points.get(order.get(i));
			Double dm = statsWithDM.points.get(order.get(i)).distancesToModel.get(dmIndex);
			accuracy.addValue(ps.getRealValue(considerPredicates), ps.predicted);
			
			if (accuracy.count > 10 && !prevDm.equals(dm))
			{
				prevDm = dm;
				adConfiguration.intervals.add(dm);
				adConfiguration.percents.add(1d * Math.round(1d * i / order.size() * 10000) / 100);
				adConfiguration.epIds.add(ps.id);
				adConfiguration.errors.add(accuracy.getAverageValue());
			}
		}
		
		return adConfiguration;
	}	
}
