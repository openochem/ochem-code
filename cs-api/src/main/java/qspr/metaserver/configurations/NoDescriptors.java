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

/**
 * Methods can be used without descriptors
 * @author itetko
 *
 */

public interface NoDescriptors {

	public void setAugmentations(Integer training, Integer validation, boolean balanced);

	public int getAugementTraining();

	public int getAugmentApply();

	public String augemenationString();

	public boolean getBalanceData();
	
	boolean isSupportAugmentation();
	
	boolean isInternalAugmentation();

	boolean isSanitize();

	void setSanitize();
	
	/**
	 * 
	 * @param smiles
	 * @param training
	 * @return null if SMILES is good
	 */
	public String isGoodSMILES(String smiles, boolean training);

	public String isPeptide(String augSmile); //added to work with peptides; to be deleted if not required

}
