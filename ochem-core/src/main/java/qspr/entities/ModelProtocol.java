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

package qspr.entities;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElements;
import javax.xml.bind.annotation.XmlRootElement;

import qspr.metaserver.configurations.BaggingConfiguration;
import qspr.metaserver.configurations.CrossValidationConfiguration;
import qspr.metaserver.configurations.ValidationConfiguration;

@XmlRootElement(name = "model-protocol")
public class ModelProtocol {

	/**
	 * If baggingConfiguration is null
	 * we do not perform any validation
	 */

	@XmlElements({
		@XmlElement(name = "baggingConfiguration", type = BaggingConfiguration.class),
		@XmlElement(name = "crossValidationConfiguration", type = CrossValidationConfiguration.class)})
	public ValidationConfiguration validationConfiguration = null;

	public String toString() {
		String res = "";

		if(validationConfiguration !=null){
			res= validationConfiguration instanceof BaggingConfiguration
					?"Bagging with " + validationConfiguration.ensembleSize + " sets"
							:"" + validationConfiguration.ensembleSize + "-fold cross-validation";

			switch (validationConfiguration.getBaggingType()) {
			case BaggingConfiguration.STRATIFIED:
				res += " (stratified)";
				break;
			default:
				break;
			}

			res += validationConfiguration.getInformativeName();

		}else
			res = "No validation";

		return res;
	}
}
