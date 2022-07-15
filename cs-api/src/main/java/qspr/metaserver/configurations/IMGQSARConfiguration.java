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

@XmlRootElement(name = "img-configuration")
public class IMGQSARConfiguration extends StandardNoDescriptorsConfiguration{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public Integer batch_size = 32;
	public Double learning_rate = 0.0005;
	public String backbone;
	//public Integer n_best;

	public IMGQSARConfiguration(){
		nepochs = 100;
	}

	@Override
	public void setAugmentations(Integer training, Integer validation, boolean balance) {
		augmentation = training; if(augmentation == null || augmentation == 1 || augmentation == 0) augmentation = null;
		augmentApplySet = validation; if(augmentApplySet == null || augmentApplySet == 1 || augmentApplySet == 0) augmentApplySet = null; 
		this.balance = balance?true:null;
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.IMGQSAR;
	}

	@Override
	public boolean isInternalAugmentation() {
		return true;
	}

	public String toString() {
		return  " batch=" + batch_size + " lr=" + learning_rate + 
				super.toString();
	}

	@Override
	public void setIterations(int iterations) {
		nepochs = iterations;
	}

	@Override
	public boolean isSupportAugmentation(){
		return true;
	}

	@Override
	public boolean isSupportConditions() {
		return true;
	}

	@Override
	public boolean isSupportRegression(){
		return true;
	}
}
