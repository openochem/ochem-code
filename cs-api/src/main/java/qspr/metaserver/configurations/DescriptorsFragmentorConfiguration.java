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

@XmlRootElement(name="fragmentor-configuration")
public class DescriptorsFragmentorConfiguration extends DescriptorsAbstractConfiguration
{
	private static final long serialVersionUID = 1L;

	public Integer fragmentType = 3;
	public Integer minFragmentLength = 2;
	public Integer maxFragmentLength = 4;

	public DescriptorsFragmentorConfiguration(){
	}

	public String toString()
	{
		return "length: " + this.minFragmentLength + "-" + this.maxFragmentLength; 
	}

	@Override
	public boolean requires3D() {
		return false;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.FRAGMENTS;
	}

	@Override
	DescriptorsFragmentorConfiguration setConfiguration(HttpServletRequest request) {

		minFragmentLength = Integer.valueOf(request.getParameter("min-length"));
		maxFragmentLength = Integer.valueOf(request.getParameter("max-length"));
		fragmentType = Integer.valueOf(request.getParameter("fragment-type"));
		return this;
	}
}
