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
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.entities.Mapping2;

@Entity
@XmlRootElement(name = "mmpair")
public class MMPair
{
	@Id
	@Column(name = "mmp_id")
	@GeneratedValue
	public Long id;
	
	@ManyToOne
	@JoinColumn(name="mol1")
	public Mapping2 molecule1;
	
	@ManyToOne
	@JoinColumn(name="mol2")
	public Mapping2 molecule2;
	
	@ManyToOne
	@JoinColumn(name="transformation_id")
	@XmlTransient
	public MMPTransformation transformation;
	
	@Column
	public Short similarity;
	
	@Transient
	public String[] values;
	
	@XmlElement
	private MMPTransformation getTransformation() {
		if (Globals.getMarshallingOption(MarshallingOption.PAIR_TRANSFORMATION))
			return transformation;
		return null;
	}
	
	public static MMPair get(Long id) {
		return (MMPair) Globals.session().get(MMPair.class, id);
	}
	
	public static MMPair thin(Long id, Integer mol1, Integer mol2) {
		MMPair mp = new MMPair();
		mp.id = id;
		mp.molecule1 = new Mapping2();
		mp.molecule2 = new Mapping2();
		mp.molecule1.id = mol1;
		mp.molecule2.id = mol2;
		
		return mp;
	}
	
	public String toString(){
		return "id="+id+ " mol1="+molecule1.id+ " mol2="+molecule2.id+ " transformation="+transformation+ " similarity="+similarity;
	}
	
}
