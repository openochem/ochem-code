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

abstract public class StandardNoDescriptorsConfiguration extends NoDescriptorsConfiguration implements EarlyStopping, Stereochemistry {

	public Double early = 0.1;
	public Boolean chirality = true;
	public Boolean shuffle;
	public Integer nepochs = 25;

	public enum ACTIVATE{
		RELU, CRELU, LRELU, ELU, SWISH, TANHEXP
	}; 

	public ACTIVATE activationFunction;

	private static final long serialVersionUID = 1L;

	@Override
	boolean areDefault() {
		return early == 0.1 && chirality == true && shuffle == null && nepochs == 25 && super.areDefault() && activationFunction == null;
	}

	@Override
	public double getEarlyStoppingFraction() {
		return early;
	}

	public boolean requiresStereochemistry(){
		return chirality == null?true:chirality;
	}

	@Override
	public boolean shuffleAllData() {
		return shuffle != null && shuffle;
	}

	@Override
	public void setIterations(int iterations) {
		nepochs = iterations;
	}

	public String toString() {
		return  (requiresStereochemistry()? " 3D":"")  + " epochs=" + nepochs +  (shuffleAllData()?" shuffle=":"") + " early=" + early +
				(activationFunction != null ?" activation=" + getActivationFunction():"") +
				super.toString();
	}

	@Override
	public String getInformativeName() {
		return super.getInformativeName() + (isSupportAugmentation() && !isInternalAugmentation()?(shuffleAllData()?"(T)":"(F)"):"") + 
				(requiresStereochemistry()?" (3D)":"");	
	}

	public String getActivationFunction() {
		String s = "" + (activationFunction ==null? ACTIVATE.RELU:activationFunction);
		return s.toLowerCase();
	}

}
