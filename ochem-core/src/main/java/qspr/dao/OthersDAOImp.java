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

package qspr.dao;

import java.util.List;

import qspr.business.toxalert.AlertsFilter;
import qspr.business.toxalert.ScreeningProcessor;
import qspr.entities.SubstructureAlert;
import qspr.metaserver.configurations.DescriptorsStructuralAlertsConfiguration;

public class OthersDAOImp implements OthersDAO {

	@Override
	@SuppressWarnings("unchecked")
	public DescriptorsStructuralAlertsConfiguration getEFGConfiguration() throws Exception {
		AlertsFilter filter = new AlertsFilter();
		filter.endpointId = 949l;
		filter.order = "id";

		List<SubstructureAlert> alerts = filter.filterCriteria().list();

		DescriptorsStructuralAlertsConfiguration conf = new DescriptorsStructuralAlertsConfiguration();
		conf.compactMode = false;
		conf.alertPatterns = ScreeningProcessor.getSMARTs(alerts);

		return conf;
	}

}
