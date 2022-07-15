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

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

@Entity
@XmlRootElement(name = "descriptor")
public class CalculatedDescriptor 
{
	@Id
	@GeneratedValue
	@Column(name = "calculated_descriptor_id")
	@XmlAttribute
	public Long id;
	
	@Column
	@XmlAttribute
	public String name;
	
	@ManyToOne
	@JoinColumn(name = "session_id")
	@XmlTransient
	public Session session;
	
	@Column
	public Double mean;
	
	@Column
	@XmlAttribute
	public Boolean selected = false;

	@Column
	@XmlAttribute
	public CalculatedDescriptorStatus status = CalculatedDescriptorStatus.NORMAL;

	public CalculatedDescriptor()
	{
		
	}
	
	public enum CalculatedDescriptorStatus {NORMAL, EXPECTED_FOUND, EXPECTED_NOT_FOUND} 
	
}
