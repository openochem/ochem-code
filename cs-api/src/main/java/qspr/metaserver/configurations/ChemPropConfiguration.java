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

@XmlRootElement(name = "chemprop-configuration")
public class ChemPropConfiguration extends NoDescriptorsConfiguration{

	private static final long serialVersionUID = 1L;

	public Integer depth = 3; // number of message passing steps
	public Integer nepochs = 100; // number of message passing steps
	public Integer hidden = 300; // hidden size of the neural network layers.
	public Integer batch = 50; //batch size

	@Override
	public String toString(){
		return "depth: " + depth + " epoch: " + nepochs + " hidden=" + hidden + //(edges != null && edges? " edges":"") + 
				" batch =" + batch + super.toString();
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.CHEMPROP;
	}

	@Override
	public void setIterations(int iterations) {
		nepochs = iterations;
	}

	public boolean isSupportRegression(){
		return true;
	}
}
