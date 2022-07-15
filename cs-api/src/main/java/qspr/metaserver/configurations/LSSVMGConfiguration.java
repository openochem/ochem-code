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

@XmlRootElement(name = "lssmv-configuration")
public class LSSVMGConfiguration extends MultiLearningAbstractConfiguration{

	private static final long serialVersionUID = 1L;

	public String kernel = "rbf";
	public Integer cv = 5;
	public Boolean useGlobaL;
	public Integer iterations;

	public String additionalParam;

	public String toString(){
		return "kernel:" + kernel + " cv:" + cv + (isGlobal() ? " global" :"" + (additionalParam == null ?"": " " + 
				((iterations != null)?" iterations = " + iterations +" ": "") + 
				additionalParam )) + super.toString();
	}

	public boolean isGlobal() {
		return useGlobaL != null && useGlobaL;
	}

	@Override
	public boolean isSupportRegression() {
		return true;
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.LSSVMG;
	}

	@Override
	public void setIterations(int iterations) {
		this.iterations = iterations;	
	}

	@Override
	public boolean isLarge(){
		return true;
	}

}
