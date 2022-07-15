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

package qspr.export;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.servlet.http.HttpServletRequest;

/**
 * A configuration fort ExportableSet class
 * @author novserj
 */

public class ExportableSetConfiguration
{
	//	private static transient final Logger logger = Logger.getLogger(ExportableSetConfiguration.class);
	public boolean useShuffleKey = false;
	public String shuffleKey;
	public boolean selectAll = false;
	public boolean useDefaults = false;
	public Map<Long, Long> propertyToUnitMap = new HashMap<Long, Long>();
	public Set<String> potentialColumnNames = new HashSet<String>();

	//TODO: Move somewhere appropriate
	@SuppressWarnings("unchecked")
	public static ExportableSetConfiguration configureFromDialog(HttpServletRequest req)
	{
		ExportableSetConfiguration c = new ExportableSetConfiguration();
		c.potentialColumnNames.addAll(req.getParameterMap().keySet());

		c.useShuffleKey = (req.getParameter("use-shuffle-key") != null);
		c.shuffleKey = req.getParameter("shuffle-key");

		for (Object param : req.getParameterMap().keySet())
		{
			String paramName = (String) param;
			if (paramName.startsWith("unit-"))
			{
				Long propertyId = Long.valueOf(paramName.substring(5));
				Long unitId = Long.valueOf(req.getParameter(paramName));
				c.propertyToUnitMap.put(propertyId, unitId);
			}
		}
		return c;
	}
}