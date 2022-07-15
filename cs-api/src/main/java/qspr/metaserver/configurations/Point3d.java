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

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;


@XmlRootElement
public class Point3d implements Serializable 
{
	private static final long serialVersionUID = 1L;
	
	@XmlTransient
	public double x, y, z;
	
	@XmlElement
	public double getX() {
		return x;
	}

	public void setX(double x) {
		this.x = x;
	}
	@XmlElement
	public double getY() {
		return y;
	}

	public void setY(double y) {
		this.y = y;
	}
	@XmlElement
	public double getZ() {
		return z;
	}

	public void setZ(double z) {
		this.z = z;
	}

	public Point3d(double x, double y, double z)
	{
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public Point3d()
	{
		this.x = 0;
		this.y = 0;
		this.z = 0;
	}
	
	public double distance(Point3d p)
	{
		return Math.sqrt((p.x-x)*(p.x-x) + (p.y-y)*(p.y-y) + (p.z-z)*(p.z-z));
	}
	
	public String toString(){
		return "" + x + ";" + y + ";" +z;
	}
}
