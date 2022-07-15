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

@XmlRootElement(name = "gin-configuration")
public class GNNConfiguration extends NoDescriptorsConfiguration{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;

	public Integer batch = 32;
	public Integer nepochs = 200;
	public Integer dim = 95;
	public Double early;

	public GNN gnn = GNN.GIN;

	public Integer patienceearly = 40;
	public Integer patience = 10;
	public Double lr_decay = 0.5;

	public enum GNN {GIN, GAIN, GGRNet};

	@Override
	public String getDefaultName() {
		return QSPRConstants.GNN;
	}

	@Override
	public String getInformativeName() {
		return super.getInformativeName() + " " + (gnn != null?gnn: "");
	}


	public boolean isSupportRegression(){
		return true;
	}

	public String toString() {
		return gnn + " batch=" + batch + " epochs=" + nepochs + (early!=null && early >0 ?" early=" + early + " ":"") + 
				" patience=" + patience + " early patience=" + patienceearly + " LR decay=" + lr_decay + " dim=" + dim +
				super.toString();
	}

	@Override
	public void setIterations(int iterations) {
		nepochs = iterations;
	}

}
