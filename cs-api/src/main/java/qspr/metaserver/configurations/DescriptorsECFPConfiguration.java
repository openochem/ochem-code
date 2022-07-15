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

@XmlRootElement(name = "ecfp-configuration")
public class DescriptorsECFPConfiguration extends DescriptorsRDKITConfiguration
{
	private static final long serialVersionUID = 1L;

	public DescriptorsECFPConfiguration(){
		super(512);
	}
	
	@Override
	public boolean requires3D() {
		return false;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.ECFP;
	}
	
	@Override
	DescriptorsAbstractConfiguration setConfiguration(HttpServletRequest request) {
		MORGAN_NBITS = request.getParameter("ecfp_nbits")!=null? Integer.valueOf(request.getParameter("ecfp_nbits")):1024;
		MORGAN_RADIUS = request.getParameter("ecfp_radius")!=null? Integer.valueOf(request.getParameter("ecfp_radius")):2;
		MORGAN_FCFP = request.getParameter("ecfp_fcfp")!=null ? true : null;
		MORGAN_COUNTS = request.getParameter("ecfp_counts")!=null ? true : false;
		dragonBlocks = 1<<9;
		return this;
	}
	
}
