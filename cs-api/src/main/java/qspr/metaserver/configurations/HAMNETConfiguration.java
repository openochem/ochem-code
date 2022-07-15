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

@XmlRootElement(name = "hamnet-configuration")
public class HAMNETConfiguration extends StandardNoDescriptorsConfiguration{

	private static final long serialVersionUID = 1L;

	/*
	public Integer f_dim = 32;
	public Integer c_dims[] = {200,200};
	public Integer he_dims = 200;
	public Integer m_radius = 4;
	public Integer mlp_dim = 100;
	*/
	
	public Boolean refine = false;
	public Integer batch_size = 64;
	public Double  learningRate = 0.001;
	public Double  dropout = 0.3;
	
	@Override
	public String getDefaultName() {
		return QSPRConstants.HAMNET;
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
	public String toString() {
		return  " learning_rate= " + learningRate +   (isRefined()?" Refine":"")
				 +  super.toString();
	}

	@Override
	public boolean isSupportAugmentation(){
		return false;
	}

	public boolean isRefined() {
		return refine != null && refine;
	}

	@Override
	public double getEarlyStoppingFraction() {
		return early;
	}


}