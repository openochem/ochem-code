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

package qspr.frontend;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

import qspr.entities.Basket;
import qspr.modelling.SetStatistics;

@XmlRootElement
public class ModelProfileData
{
	public List<SetStatsProfile> setStats = new ArrayList<SetStatsProfile>();
	
	public void addStats(Basket b, SetStatistics stats, SetStatistics estimatedStats) {
		SetStatsProfile ssProfile = new SetStatsProfile();
		ssProfile.basket = b;
		ssProfile.statistics = stats;
		ssProfile.estimatedStatistics = estimatedStats;
		setStats.add(ssProfile);
	}
}

class SetStatsProfile {
	public Basket basket;
	public SetStatistics statistics;
	public SetStatistics estimatedStatistics;
	
}
