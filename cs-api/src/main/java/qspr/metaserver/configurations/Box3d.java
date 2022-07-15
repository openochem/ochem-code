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
public class Box3d implements Serializable 
{
	private static final long serialVersionUID = 1L;
	
	@XmlTransient
	public Point3d center, size;
	
	public Box3d(Point3d center, Point3d size)
	{
		this.center = center;
		this.size = size;
	}
	@XmlElement
	public Point3d getCenter() {
		return center;
	}

	public void setCenter(Point3d center) {
		this.center = center;
	}
	@XmlElement
	public Point3d getRowsSize() {
		return size;
	}

	public void setSize(Point3d size) {
		this.size = size;
	}

	public Box3d()
	{
	}

	public String toString(){
		return center + "/" + size; 
	}
}
