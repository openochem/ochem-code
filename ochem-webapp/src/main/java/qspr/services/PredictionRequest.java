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

package qspr.services;

public class PredictionRequest {
	public long modelId;
	
	public String[] sdfs;
	
	public String predictionScenario;
	
	public String sessionGUID;
	
	public long basketId;

	public long getModelId() {
		return modelId;
	}

	public void setModelId(long modelId) {
		this.modelId = modelId;
	}

	public String[] getSdfs() {
		return sdfs;
	}

	public void setSdfs(String[] sdfs) {
		this.sdfs = sdfs;
	}

	public String getPredictionScenario() {
		return predictionScenario;
	}

	public void setPredictionScenario(String predictionScenario) {
		this.predictionScenario = predictionScenario;
	}

	public String getSessionGUID() {
		return sessionGUID;
	}

	public void setSessionGUID(String sessionGUID) {
		this.sessionGUID = sessionGUID;
	}

	public long getBasketId() {
		return basketId;
	}

	public void setBasketId(long basketId) {
		this.basketId = basketId;
	}
}
