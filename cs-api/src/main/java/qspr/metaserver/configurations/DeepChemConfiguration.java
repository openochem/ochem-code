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

import java.util.HashSet;
import java.util.Set;

import javax.xml.bind.annotation.XmlRootElement;

import qspr.workflow.utils.QSPRConstants;

@XmlRootElement(name = "deepchem-configuration")
public class DeepChemConfiguration extends StandardNoDescriptorsConfiguration{
	private static final long serialVersionUID = 1L;

	public DeepChemMethod method;

	public enum DeepChemMethod
	{
		DAG, GRAPH_CONV, MPNN,  TEXTCNN, WEAVE;

		public static boolean isRegressionMethod(DeepChemMethod method){
			return method == null || method == GRAPH_CONV || method == DAG || method == MPNN 
					|| method == TEXTCNN || method == WEAVE //|| method == MULTITASK 
					;
		}

	}

	/*public Double learning_rate = 0.001;
	public Double dropout = 0.25;
	public Integer dense_layer_size = 128;
	public String  graph_conv_layers = "64,64";
	public Integer M = 5;
	public Integer T = 3;
	public Integer n_hidden = 100;
	public Integer n_embedding = 75;	
	 */

	public Integer maxlength; // max length of the training set sequences
	Set<Character> chars = new HashSet<Character>();

	@Override
	public boolean isSupportRegression(){
		return true;
	}

	@Override
	public String toString(){
		return getDefaultName() + super.toString();
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.DEEPCHEM;
	}

	@Override
	public String getInformativeName() {
		return super.getInformativeName() +  (method == null ? "" :  method);	
	}

	public boolean isTEXTCNN() {
		return method == DeepChemMethod.TEXTCNN;
	}

	@Override
	public boolean isSupportAugmentation(){
		return isTEXTCNN();
	}

	@Override
	public void setIterations(int iterations) {
		nepochs = iterations;
	}

	@Override
	public String isGoodSMILES(String smiles, boolean training) {
		if(isTEXTCNN()) {
			if(training) {
				if(maxlength == null || smiles.length()>maxlength)maxlength = smiles.length();
				for (int i = 0;i < smiles.length(); i++)
					chars.add(smiles.charAt(i));
			}else {
				if(maxlength != null && smiles.length()>maxlength*1.1)return " too large SMILES > " + maxlength;
				if(chars != null && chars.size() > 0) // compatibility with previous codes
					for (int i = 0;i < smiles.length(); i++)
						if(!chars.contains(smiles.charAt(i)))return " contains non - supported character " + smiles.charAt(i);
			}
		}
		return super.isGoodSMILES(smiles, training);
	}

	@Override
	public boolean isSanitize() {
		return sanitize != null && sanitize && isTEXTCNN();
	}

}
