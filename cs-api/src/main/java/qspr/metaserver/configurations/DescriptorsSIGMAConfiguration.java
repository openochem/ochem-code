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

import com.eadmet.exceptions.UserFriendlyException;

@XmlRootElement(name = "sigma-configuration")
public class DescriptorsSIGMAConfiguration extends DescriptorsAbstractConfiguration{
	private static final long serialVersionUID = 1L;

	public Double sigmawidth;
	final static double  DEFAULT_SIGMA = 0.001;
	final static public double  OPTION_WIDTHS []= {0.005, 0.001, 0.0005};
	public Boolean all;

	@Override
	public boolean requires3D() {
		return true;
	}

	@Override
	public boolean isCachable() {
		return true;
	}

	@Override
	public boolean isLongCalculation()
	{
		return true;
	}

	@Override
	public String getDefaultTypeName() {
		return DescriptorsConfiguration.SIGMA;
	}

	@Override
	public DescriptorsAbstractConfiguration setAllOn(){
		all = true;
		sigmawidth = null;
		return this;
	}

	public double getWidth() {
		return sigmawidth == null?DEFAULT_SIGMA:sigmawidth;
	}

	@Override
	public String toString() {
		return "width="+getWidth();
	}

	@Override
	DescriptorsAbstractConfiguration setConfiguration(HttpServletRequest request) {
		sigmawidth = Double.valueOf(request.getParameter("sigmawidth"));

		boolean found = false;
		for(double v:OPTION_WIDTHS)
			if(v == sigmawidth)found = true;

		if(!found)throw new UserFriendlyException("sigma bin width is not in the allowed range");

		if(sigmawidth == DEFAULT_SIGMA)
			sigmawidth = null;
		return this;
	}

}
