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

import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "oestate-configuration")
public class DescriptorsOEstateConfiguration extends DescriptorsAbstractConfiguration
{
	private static final long serialVersionUID = 2L;

	public int count = 0;
	public int bond = 1;
	public String info;

	public String toString()
	{
		return "";
	}

	@Override
	public boolean requires3D() {
		return false;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.OEstate;
	}

	@Override
	public DescriptorsAbstractConfiguration setAllOn() {
		bond = 1;
		return this;
	}

	@Override
	DescriptorsAbstractConfiguration setConfiguration(HttpServletRequest request) {
		count = request.getParameter("oestate_count") != null ? 1 : 0;
		bond = request.getParameter("oestate_bond") != null ? 1 : 0;
		return this;
	}
}
