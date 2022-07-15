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

@XmlRootElement(name = "attfp-configuration")
public class ATTFPConfiguration extends NoDescriptorsConfiguration implements SupportsOneOutputOnly, MaximalSizeRestriction{

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;


	public Integer batch = 64;
	public Integer nepochs = 200;
	public Integer patience_reduce = 2;
	public Integer patience_early = 40;
	public Integer radius = 2;
	public Integer T = 2;
	public Integer fp_dim = 200;
	public Integer cosineT = 14;

	public Double lr = 0.0032;
	public Double dropout = 0.2;
	public Double weight_decay = 0.00001;

	public Boolean cosine = true;
	public Boolean early = true;

	public Boolean simpleO;
	public Boolean singleT;
	public Boolean lngru;

	public Integer molsize;

	@Override
	public int getMaxSize() {
		return molsize == null? 200: molsize;
	}
	
	@Override
	public String getDefaultName() {
		return QSPRConstants.ATTFP;
	}

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

}
