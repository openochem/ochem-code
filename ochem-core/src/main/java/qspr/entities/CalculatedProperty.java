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

import java.text.DecimalFormat;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlValue;

import qspr.Globals;

@Entity
public class CalculatedProperty 
{
	@Id
	@GeneratedValue
	@Column(name = "calc_property_id")
	@XmlTransient
	private Long id;
	
	@XmlTransient
	@ManyToOne
	@JoinColumn(name = "molecule_id")
	public Molecule molecule;
	
	@ManyToOne
	@JoinColumn(name = "property_id")
	@XmlTransient
	public Property property;
	
	@Column
	@XmlAttribute
	public Double accuracy;
	
	@XmlAttribute@Transient
	public String model;
	
	@XmlAttribute@Transient
	public String error;
	
	@Transient@XmlTransient
	private String name;
	
	@Transient@XmlTransient
	public boolean dontDelete = false;
	
	static public CalculatedProperty get(Molecule molecule, Property property)
	{
		return null;
	}
	
	@XmlAttribute
	public String getName()
	{
		if (property != null)
			return property.getName();
		else
			return name;
	}
	
	public void setName(String _name)
	{
		name = _name;
		if (property == null)
			property = new Property();
		property.setName(_name);
	}
	
	public void setModel(String _model)
	{
		model = _model;
	}
	
	@XmlValue
	public double value;
	
	@XmlAttribute(name = "printable-value")
	public String getPrintableValue()
	{
		if (value == Globals.ERROR_RESULT)
			return "Error";
		return (new DecimalFormat("#0.000")).format(value);
	}
	
	@XmlAttribute(name = "id")
	public Long getPropertyId()
	{
		if (property != null)
			return property.id;
		else
			return null;
	}
	
	public void setPropertyId(Long _id)
	{
		if (property == null)
			property = new Property();
		property.id = _id;
	}
}
