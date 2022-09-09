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

import javax.xml.bind.annotation.XmlRootElement;

import qspr.workflow.utils.QSPRConstants;

@XmlRootElement(name = "dimenet-configuration")
public class DIMENETConfiguration extends NoDescriptorsConfiguration implements Supports3D{

	private static final long serialVersionUID = 1L;

	public Integer batch = 32;
	public Integer nbepochs = 100;
	public Boolean early = true;
	public Boolean external3D;

	@Override
	public String getDefaultName() {
		return QSPRConstants.DIMENET;
	}


	@Override
	public boolean isSupportDescriptors() {
		return false;
	}

	@Override
	public boolean isSupportRegression(){
		return true;
	}


	@Override
	public String toString(){
		return super.toString() + 
				(" batch=" +  batch) + (" epochs=" + nbepochs);  
	}

	@Override
	public void setIterations(int iterations) {
		nbepochs = iterations;
	}

	@Override
	public boolean isSupportAugmentation(){
		return false;
	}

	@Override
	public boolean requires3D() {
		return external3D != null && external3D;
	}

	@Override
	public void setUse3D(boolean yes) {
		if(yes)external3D = true;
		else
			external3D = null;
	}
}