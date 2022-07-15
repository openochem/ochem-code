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

package qspr.business.toxalert;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import qspr.Globals;
import qspr.entities.SubstructureAlert;

/**
 * An in-memory cache of all substructure alerts used to speed up page loading
 * @author midnighter
 *
 */
public class AlertPool
{
	private List<SubstructureAlert> allAlerts = null;
	private Map<String, SubstructureAlert> alertsByURL = new HashMap<String, SubstructureAlert>();

	/**
	 * Intended for caching the alerts list to increase page responsiveness
	 * @return
	 */
	@SuppressWarnings("unchecked")
	public List<SubstructureAlert> getAllAlerts()
	{
		if (allAlerts != null)
			return allAlerts;
		else
		{
			allAlerts = Globals.session().createCriteria(SubstructureAlert.class).list();
			for (SubstructureAlert alert : allAlerts)
			{
				alert.urlFriendlyName = urlFriendly(alert.name);
				alertsByURL.put(alert.urlFriendlyName, alert);
				Globals.session().evict(alert);
				alert.article = null;
				alert.property = null;
			}
			return allAlerts;
		}
	}

	public SubstructureAlert getAlertByURL(String urlFriendlyName)
	{
		getAllAlerts();
		return alertsByURL.get(urlFriendlyName);
	}

	private String urlFriendly(String str)
	{
		return str.replaceAll("[^a-zA-Z\\-]+", "_").toLowerCase();
	}
}
