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

package qspr.modelling.applier;

public class Prediction
{
	private Long depictionID;
	private Integer moleculeID;
	private String error;
	private String Smiles;
	private String InChIKey;
	private PropertyPrediction[] predictions;
	
	public Long getDepictionID()
	{
		return depictionID;
	}

	public void setDepictionID(Long depictionID)
	{
		this.depictionID = depictionID;
	}

	public Integer getMoleculeID()
	{
		return moleculeID;
	}

	public void setMoleculeID(Integer moleculeID)
	{
		this.moleculeID = moleculeID;
	}
	
	public String getError()
	{
		return error;
	}
	
	public void setError(String error)
	{
		this.error = error;
	}

	public void setPredictions(PropertyPrediction[] predictions)
	{
		this.predictions = predictions;
	}

	public PropertyPrediction[] getPredictions()
	{
		return predictions;
	}

	public String getSmiles() {
		return Smiles;
	}

	public void setSmiles(String inputSmiles) {
		Smiles = inputSmiles;
	}

	public String getInChIKey() {
		return InChIKey;
	}

	public void setInChIKey(String inChIKey) {
		InChIKey = inChIKey;
	}
	
	
}
