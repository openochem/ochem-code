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
import java.io.Serializable;

import javax.persistence.Column;
import javax.persistence.Embeddable;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.IdClass;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

@Entity
@XmlRootElement(name = "moleculename")
@IdClass(ColorNamePk.class)
public class ColoredName extends AbstractMoleculeName
{
	@Id
	Long exp_property_id;
	@Id
	@XmlAttribute(name="id")
	Long molecule_name_id;
	
	
	@XmlAttribute
	public String name;
	
	/* validation encoding
	 * 1 = green
	 * 2 = dark red (considering stereochemistry
	 * 3 = blue
	 * 4 = red
	 * else black
	 */
	@Id
	@XmlAttribute
	public Long validation;
	
	@Column(name="userName")
	@XmlAttribute
	public String user;

	@ManyToOne
	@JoinColumn(name = "exp_property_id", insertable=false, updatable=false)
	@XmlTransient
	public ExperimentalProperty ep;
	
	public String toString()
	{
		return name;
	}

	@Override
	public String getName() 
	{
		return name;
	}

	@Override
	public int getValidation() 
	{
		return validation.intValue();
	}
	
	
	
}

@Embeddable
class ColorNamePk implements Serializable
{
	private static final long serialVersionUID = 1L;
	Long exp_property_id;
	Long molecule_name_id;
	Long validation;
	
	public boolean equals(ColorNamePk o)
	{
		return false;
	}
	
}