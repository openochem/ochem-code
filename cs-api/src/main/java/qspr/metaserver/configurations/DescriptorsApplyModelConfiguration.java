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

import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.workflow.utils.QSPRConstants;

@XmlRootElement(name = "applymodel-configuration")
public class DescriptorsApplyModelConfiguration extends DescriptorsAbstractConfiguration {

	private static final long serialVersionUID = 1L;

	public long modelId;
	public PredictionScenario scenario;
	public String units[];

	public DescriptorsApplyModelConfiguration(){
	}

	public DescriptorsApplyModelConfiguration(long id){
		modelId = id;
	}

	@Override
	public String toString(){
		return "model:" + modelId;
	}

	@Override
	public boolean requires3D() {
		return false;
	}

	@Override
	public boolean isCachable(){
		return modelId < QSPRConstants.MAXIMAL_PUBLIC_ID; 
	}

	@Override
	public void cleanXML(){
		super.cleanXML();
		if(scenario == PredictionScenario.PREDICTION_ONLY)scenario = null;
	}

	public boolean keepColumnsInOrder() {
		return true;
	}

	@Override
	public String getDefaultTypeName(){
		return DescriptorsConfiguration.ApplyModel + "_" + modelId;
	}

}
