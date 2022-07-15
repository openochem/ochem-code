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

@XmlRootElement(name = "epa-configuration")
public class DescriptorsEPAConfiguration extends DescriptorsAbstractConfiguration {

	private static final long serialVersionUID = 1L;

	public Boolean epa3D; // only for the future so far..

	public DescriptorsEPAConfiguration() {
		moleculeTimeout = 5;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.EPA;
	}

	@Override
	public boolean isLongCalculation(){
		return true;
	}

	@Override
	public String toString(){
		return "T.E.S.T.";
	}

	@Override
	public boolean isCachable() {
		return true;
	}

	@Override
	public DescriptorsAbstractConfiguration setAllOn() {
		epa3D = true;
		return this;
	}

	@Override
	public boolean requires3D() {
		return epa3D != null? epa3D:false;
	}

	@Override
	DescriptorsAbstractConfiguration setConfiguration(HttpServletRequest request) {
		epa3D = request.getParameter("epa3D")!=null? true : null;
		return this;
	}
}
