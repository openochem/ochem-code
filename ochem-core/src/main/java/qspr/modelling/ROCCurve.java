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

package qspr.modelling;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;



@XmlRootElement
public class ROCCurve
{
	@XmlAttribute
	public String setId;
	
	@XmlElement
	public List<ROCPoint> points = new ArrayList<ROCPoint>();
	
	public void addPoint(Double x, Double y, Long id)
	{
		points.add(new ROCPoint(x, y, id));
	}
	
	public void updatePoint(Double x, Double y, Long id)
	{
		ROCPoint point = points.get(points.size() - 1);
		point.x = x;
		point.y = y;
		point.id = id;
	}
	
	public ROCCurve()
	{
		
	}
	
	public ROCCurve(String setId)
	{
		this.setId = setId;
	}
}

@XmlRootElement
class ROCPoint
{
	@XmlAttribute
	Double x;
	@XmlAttribute
	Double y;
	@XmlAttribute
	Long id;
	
	public ROCPoint(Double x, Double y, Long id)
	{
		this.x = x;
		this.y = y;
		this.id = id;
	}
	
	public ROCPoint()
	{
		
	}
}