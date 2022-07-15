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

import java.io.IOException;

import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "mordred-configuration")
public class DescriptorsMORDREDConfiguration extends DescriptorsAbstractConfiguration{

	private static final long serialVersionUID = 1L;

	public Descr mordred = Descr.D3;

	public enum Descr {D2, D3};

	@Override
	public String toString() {
		return mordred == Descr.D2? " 2D": " All";
	}

	@Override
	protected void setAllOn2D() {
		mordred = Descr.D2;
	}

	@Override
	public boolean isLongCalculation(){
		return false;
	}

	@Override
	public boolean requires3D() {
		return mordred == Descr.D3;
	}

	@Override
	public DescriptorsAbstractConfiguration setAllOn() throws IOException{
		mordred=Descr.D3;
		return this;
	}

	@Override
	public boolean isCachable() {
		return false;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.MORDRED;
	}

	@Override
	DescriptorsMORDREDConfiguration setConfiguration(HttpServletRequest request) {
		if (request.getParameter("mordred") != null) {
			if(Descr.valueOf(request.getParameter("mordred")) == Descr.D2) mordred = Descr.D2;
			else
				mordred = Descr.D3;
		}
		return this;	
	}

}
