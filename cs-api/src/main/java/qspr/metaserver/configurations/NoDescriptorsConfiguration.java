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

import qspr.dao.Various;
import qspr.workflow.utils.QSPRConstants;

abstract public class NoDescriptorsConfiguration extends MultiLearningAbstractConfiguration implements NoDescriptors{

	private static final long serialVersionUID = 1L;

	public Integer augmentation;
	public Integer augmentApplySet;
	public Boolean balance;
	public Boolean asPeptides;
	public Boolean sanitize; // default -- do not sanitize

	public NoDescriptorsConfiguration() {
		scaleTypeY = ScalingType.STANDARDIZE;
		scaleTypeX = ScalingType.STANDARDIZE;
	}

	boolean areDefault() {
		return augmentation == null && augmentApplySet == null && balance == null && asPeptides == null;
	}

	@Override
	public void setAugmentations(Integer training, Integer validation, boolean balance) {
		augmentation = training; if(augmentation == null || augmentation <= 1) augmentation = null;
		augmentApplySet = validation; if(augmentApplySet == null || augmentApplySet <=1) augmentApplySet = null;
		this.balance = balance?true:null;
	}

	@Override
	final public int getAugementTraining() {
		return !isSupportAugmentation() || augmentation == null ? 1: augmentation;
	}

	@Override
	final public int getAugmentApply() {
		return !isSupportAugmentation() || augmentApplySet == null?1: augmentApplySet;
	}

	@Override
	public String toString(){
		return  augemenationString() + " " + (isSanitize()?" sn ":"") +  super.toString();
	}

	@Override
	public  String augemenationString() {
		String str = !isSupportAugmentation() || ( getAugmentApply() == 1  && getAugementTraining() == 1) ?"": 
			" " + getAugementTraining() + "/" + getAugmentApply();
		return balance != null && balance? str + " balanced":str;
	}

	@Override
	public boolean isSupportDescriptors() { // by default
		return false;
	}

	@Override
	public boolean getBalanceData() {
		return isSupportAugmentation()? (balance==null?false:balance):false;
	}

	@Override
	public boolean isSupportAugmentation(){
		return false;
	}

	@Override
	public boolean isInternalAugmentation() {
		return false;
	}

	@Override
	public String isGoodSMILES(String smiles, boolean training) {
		return smiles.indexOf(QSPRConstants.ERROR) == -1 ? null:smiles;
	}

	@Override
	public boolean isLarge(){
		return true;
	}

	public String isPeptide(String augSmile){
		if(asPeptides == null || !asPeptides) return null;

		try {
			return Various.molecule.convertToFormat(augSmile,QSPRConstants.PEPTIDES);
		}catch(Exception e) {
			return QSPRConstants.ERROR;
		}
	} //added to work with peptides; to be deleted if not required

	@Override
	public boolean isSanitize() {
		return sanitize != null && sanitize;
	}

	@Override
	public void setSanitize() {
		sanitize = true;
	}

	@Override
	public String getInformativeName() {
		return super.getInformativeName() + " ";	
	}
	
}
