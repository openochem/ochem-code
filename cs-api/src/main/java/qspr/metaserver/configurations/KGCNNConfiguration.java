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

import qspr.dao.Various;
import qspr.workflow.utils.QSPRConstants;

@XmlRootElement(name = "kgcnn-configuration")
public class KGCNNConfiguration extends NoDescriptorsConfiguration implements SupportsOneOutputOnly, Supports3D{

	private static final long serialVersionUID = 1L;
	public Boolean external3D;

	public enum KGCNN {

		Schnet,
		PAiNN, 
		GATv2, 
		GINE,

		GCN,
		ChemProp, 
		GraphSAGE,
		AttFP,
		GIN,
		GAT, 
		HamNet, 
		DimeNetPP, 


	};

	public KGCNN method;

	public Integer batch = 32;
	public Integer nepochs = 200;

	@Override
	public String getDefaultName() {
		return QSPRConstants.KGCNN;
	}

	public String getRealName() {
		return ""+method;	
	}

	@Override
	public String getInformativeName() {
		return super.getInformativeName() +  (method == null ? "" : getRealName());	
	}

	@Override
	public boolean isSupportRegression(){
		return true;
	}

	public String toString() {
		return  " batch=" + batch + " epochs=" + nepochs +
				super.toString();
	}

	@Override
	public void setIterations(int iterations) {
		nepochs = iterations;
	}

	@Override
	public String isGoodSMILES(String smiles, boolean training) {
		try {
			if(Various.molecule.getAtomCount(smiles)<2) return QSPRConstants.ERROR_ONE_ATOM ;
		}catch (Exception e) {}
		return super.isGoodSMILES(smiles, training);
	}

	@Override
	public boolean requires3D() {
		return external3D != null && external3D && (method == KGCNN.PAiNN || method == KGCNN.DimeNetPP || method == KGCNN.Schnet || method == KGCNN.HamNet);
	}

	@Override
	public void setUse3D(boolean yes) {
		if(yes)external3D = true;
		else
			external3D = null;
		if(!requires3D())external3D = null; // to disable for those which do not require 3D
	}

	@Override
	public boolean isSupportConditions(){
		return isSupportDescriptors();
	}

	@Override
	public boolean isSupportDescriptors() {
		return false;
	}

}
