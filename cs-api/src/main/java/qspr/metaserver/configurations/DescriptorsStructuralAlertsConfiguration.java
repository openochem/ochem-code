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

package qspr.metaserver.configurations;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

/**
 * A configuration for the StructuralAlerts tasks
 * @author midnighter
 *
 */
@XmlRootElement(name = "structuralalerts-configuration")
public class DescriptorsStructuralAlertsConfiguration extends DescriptorsAbstractConfiguration
{
	private static final long serialVersionUID = 1L;

	/**
	 * A set of alert patterns in the extended SMARTS format (see ExtendedSMART.java)
	 */
	public List<String> alertPatterns = new ArrayList<String>();

	/**
	 * Return a bit per alert: 1 for a match and 0 otherwise. Store 64 bits in 1 long value
	 */
	public boolean compactMode;

	public DescriptorsStructuralAlertsConfiguration(){
	}

	public void setAlerts(List<String> addAlerts){
		alertPatterns = addAlerts;
	}

	public String toString() { return ""; }

	public String printConfig() {

		String s = "" + alertPatterns.size() + " alert patterns are used. First several patterns are:\n";

		for (int i = 0; i < alertPatterns.size() && i < 5; i++)
			s += alertPatterns.get(i) + "; ";

		return s;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.StructuralAlerts;
	}

	public int size() {
		return alertPatterns == null?0:alertPatterns.size();
	}

	public String get(int i) {
		return alertPatterns.get(i);
	}

}
