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

public class DataReferenceRequest{

	private String dataReference;
	private String sessionGUID;
	private long modelId;
	private String predictionScenario;
	private int datasize;

	public String getSessionGUID() {
		return sessionGUID;
	}
	public void setSessionGUID(String sessionGUID) {
		this.sessionGUID = sessionGUID;
	}
	public long getModelId() {
		return modelId;
	}
	public void setModelId(long modelId) {
		this.modelId = modelId;
	}
	public String getPredictionScenario() {
		return predictionScenario;
	}
	public void setPredictionScenario(String predictionScenario) {
		this.predictionScenario = predictionScenario;
	}
	public String getDataReference() {
		return dataReference;
	}
	public void setDataReference(String dataKeyReference) {
		this.dataReference = dataKeyReference;
	}

	public String toString(){
		return "reference="+dataReference+" sessionGUID="+sessionGUID+ " modelId="+modelId+ "predictionScenario="+predictionScenario;
	}
	public int getDatasize() {
		return datasize;
	}
	public void setDatasize(int datasize) {
		this.datasize = datasize;
	}
	
}
	

