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

import qspr.metaserver.configurations.CNFConfiguration.NORMALIZE;
import qspr.workflow.utils.QSPRConstants;
@XmlRootElement(name = "eagcnг-configuration")
public class EAGCNGConfiguration extends StandardNoDescriptorsConfiguration implements MaximalSizeRestriction, Iterations{

	private static final long serialVersionUID = 1L;

	public Integer batchsize = 512;
	public Integer molsize = 100;
	public Double dropout = 0.2;
	public Double learningRate = 0.01;
	public Double weightDecay = 0.0001;
	public String method = "concate";
	public String n_sgc1 = "30,15,15,15,15";
	public String n_sgc2 = "60,30,30,30,30";
	public String n_sgc3 = null;
	public String n_den = "64,32,12";
	public Boolean gate = false;
	public NORMALIZE normalisation = NORMALIZE.LAYER;


	public EAGCNGConfiguration() {
		if(areDefault())
			nepochs = 1000;
	}

	@Override
	public boolean isSupportRegression() {
		return true;
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.EAGCNG;
	}

	@Override
	public String toString() {
		return " epochs=" + nepochs + " norm=" + normalisation +
				" gate=" + gate + " batch=" +  batchsize + " dropout=" + dropout + " learningRate=" + learningRate + " weightDecay=" + weightDecay +
				" method=" + method + " n_sgc1=" + n_sgc1 +
				((n_sgc1.equalsIgnoreCase(n_sgc2) && n_sgc1.equalsIgnoreCase(n_sgc3))
						? (" n_sgc1 = n_sgc2 = n_sgc3") 
								: ((n_sgc2 != null ? n_sgc1.equalsIgnoreCase(n_sgc2)? " n_sgc2 = n_sgc1" : " n_sgc2=" + n_sgc2 :"") + 
										(n_sgc3 != null ? n_sgc3.equalsIgnoreCase(n_sgc1)? " n_sgc3 = n_sgc1" : " n_sgc3=" + n_sgc3 :""))) 
				+ " n_den=" +  n_den;
	}

	@Override
	public int getMaxSize() {
		return molsize == null? 1000: molsize;
	}

	@Override
	public void setIterations(int iterations) {
		nepochs = iterations;
	}

}
