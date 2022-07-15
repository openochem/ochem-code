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

import javax.xml.bind.annotation.XmlRootElement;

import qspr.modelling.applier.Prediction;

@XmlRootElement
public class ModelResponse 
{
	private String status;
	private String detailedStatus;


	private Prediction[] predictions;
	private String modelDescriptionUrl;
	private long taskId;
	private long metaserverTaskId;
	
	public long getMetaserverTaskId() {
		return metaserverTaskId;
	}

	public void setMetaserverTaskId(long metaserverTaskId) {
		this.metaserverTaskId = metaserverTaskId;
	}

	public long getTaskId()
	{
		return taskId;
	}
	
	public void setTaskId(long taskId)
	{
		this.taskId = taskId;
	}
	
	public String getDetailedStatus() 
	{
		return detailedStatus;
	}

	public void setDetailedStatus(String detailedStatus) 
	{
		this.detailedStatus = detailedStatus;
	}
	
	/**
	 * @return the status
	 */
	public String getStatus() 
	{
		return status;
	}
	
	/**
	 * @param status the status to set
	 */
	public void setStatus(String status) 
	{
		this.status = status;
	}
	
	/**
	 * @return the predictionValue
	 */
	public Prediction[] getPredictions() 
	{
		return predictions;
	}
	
	/**
	 * @param predictionValue the predictionValue to set
	 */
	public void setPredictions(Prediction[] predictions) 
	{
		this.predictions = predictions;
	}
	
	/**
	 * @return the modelDescriptionUrl
	 */
	public String getModelDescriptionUrl() 
	{
		return modelDescriptionUrl;
	}
	
	/**
	 * @param modelDescriptionUrl the modelDescriptionUrl to set
	 */
	public void setModelDescriptionUrl(String modelDescriptionUrl) 
	{
		this.modelDescriptionUrl = modelDescriptionUrl;
	}
}
