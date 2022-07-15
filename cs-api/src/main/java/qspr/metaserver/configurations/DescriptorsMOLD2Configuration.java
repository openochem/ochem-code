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

@XmlRootElement(name="mold2-configuration")
public class DescriptorsMOLD2Configuration extends DescriptorsAbstractConfiguration
{
	private static final long serialVersionUID = 1L;

	public DescriptorsMOLD2Configuration()
	{

	}

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
		return DescriptorsConfiguration.MOLD2;
	}

	@Override
	DescriptorsMOLD2Configuration setConfiguration(HttpServletRequest request) {
		return this;
	}
}
