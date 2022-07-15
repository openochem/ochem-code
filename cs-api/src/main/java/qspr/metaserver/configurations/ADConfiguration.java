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

import java.io.Serializable;
import java.util.List;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "ad-configuration")
public class ADConfiguration implements Serializable
{
	private static final long serialVersionUID = 2L;
	
	//Some test comment
	public String dmName;
	public int pointsPerBlock = 20;
	public double windowSizeInPercent = 0.1;
	public boolean allowErrorFall = false;
	public String averagingType;
	
	@XmlElement(name = "interval")
	public List<Double> intervals;
	
	@XmlElement(name = "error")
	public List<Double> errors;
	
	@XmlElement
	public List<Double> percents;
	
	@XmlElement(name = "epId")
	public List<Long> epIds;
	
	public Double threshold;
	
	
	/**
	 *  Provides predicted error (RMSE) for regression model
	 *  Or unbalanced accuracy (i.e., % of correct predictions) for classification model
	 *  as function of the dm value
	 * @param dm
	 * @return
	 */
	
	public double getAverageAccuracy(double dm)
	{
		int interval = 0;
		while (interval < intervals.size() - 1 && intervals.get(interval) < dm)
			interval++;
		return errors.get(interval);
	}
	
	public ADConfiguration setAllowErrorFall(boolean allowErrorFall)
	{
		this.allowErrorFall = allowErrorFall;
		return this;
	}
	
	public ADConfiguration()
	{
		
	}
	
	public ADConfiguration(int pointsPerBlock)
	{
		this.pointsPerBlock = pointsPerBlock;
	}
	
}
