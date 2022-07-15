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

@XmlRootElement(name = "dnn-configuration")
public class DNNConfiguration extends MultiLearningAbstractConfiguration{

	private static final long serialVersionUID = 2L;

	public String modeltype = "dense7new";

	public String optimizertype = "adam";

	public Integer epochs = 1000;

	public Integer batchSize = 1024;

	public Float ratio = 0.8f;

	public String additionalParam;

	public String additionalParamModel;

	public String activation;

	public DNNConfiguration() {
		scaleTypeY = ScalingType.STANDARDIZE;
		scaleTypeX = ScalingType.STANDARDIZE;
	}

	@Override
	public String toString()
	{
		return super.toString() +  " DNN: " + modeltype + " epochs: " + epochs 
				+ " optimizer: " + optimizertype + " batch: " + batchSize 
				+ (activation != null ? " activation: "+activation:"");
	}

	@Override
	public boolean isSupportRegression() {
		return true;
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.DNN;
	}

	@Override
	public boolean isSupportPredicates() {
		return true;
	}

	@Override
	public void setIterations(int iterations) {
		epochs = iterations;
	}

	/*
	@Override
	public int requireMinimumRecords() {
		return 64;
	}
	 */

	@Override
	public boolean isLarge(){
		return true;
	}

}
