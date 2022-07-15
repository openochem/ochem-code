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

@XmlRootElement(name = "macau-configuration")
public class MACAUConfiguration extends MultiLearningAbstractConfiguration {

	private static final long serialVersionUID = 1L;

	public Double test_percent = 0.1;
	public Integer num_latent = 32;
	public Double accuracy = 0.5;
	public Integer burnin  = 512;
	public Integer samples  = 1600;
	public Boolean adaptive;

	public String toString(){
		return  "ratio=" + test_percent + " latent=" + num_latent + (isAdaptive()?" adaptive":(" accuracy=" + accuracy)) + " burnin=" + burnin + " samples=" + samples + 
				super.toString();
	}

	@Override
	public boolean isSupportRegression() {
		return true;
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.MACAU;
	}

	@Override
	public void setIterations(int iterations) {
	}

	public boolean isAdaptive() {
		return adaptive != null && adaptive;
	}

	@Override
	public boolean isLarge(){
		return true;
	}

}
