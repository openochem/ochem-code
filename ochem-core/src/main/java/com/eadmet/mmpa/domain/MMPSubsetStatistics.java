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

package com.eadmet.mmpa.domain;

import javax.persistence.Column;
import javax.persistence.Embeddable;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import com.eadmet.utils.NumericalValueStandardizer;

/**
 * Statistics over some subset of MMPs.
 * Typically, of MMPs from a particular transformation and particular dataset.
 * 
 * @author midnighter
 */
@XmlRootElement
@Embeddable
public class MMPSubsetStatistics
{
	@Column
	@XmlElement
	public double deltaMean;
	
	@Column
	@XmlElement
	public double deltaStd;
	
	@Column
	@XmlElement
	public int pairsCount;
	
	@XmlAttribute
	@Column
	public Integer nPP;
	
	@XmlAttribute
	@Column
	public Integer nNN;
	
	@XmlAttribute
	@Column
	public Integer nNP;
	
	@XmlAttribute
	@Column
	public Integer nPN;
	
	@Column
	@XmlElement
	public double pValue;
	
	public MMPSubsetStatistics() {
		
	}
	
	public void invert()
	{
		deltaMean = -deltaMean;
		Integer tmp = nNP;
		nNP = nPN;
		nPN = tmp;
	}
	
	public MMPSubsetStatistics(double deltaMean, double deltaStd, int pairsCount) {
		this.deltaMean = NumericalValueStandardizer.getSignificantDigitsDouble(deltaMean, 2);
		this.deltaStd = NumericalValueStandardizer.getSignificantDigitsDouble(deltaStd, 2);
		
		// A special case - STD for only one sample is NaN and should be zero
		if (Double.isNaN(this.deltaStd))
			this.deltaStd = 0;
		
		this.pairsCount = pairsCount;
	}
}
