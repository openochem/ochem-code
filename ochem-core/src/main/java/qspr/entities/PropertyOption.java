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
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.Globals;
import qspr.ThreadScope;

@Entity
@XmlRootElement(name = "option")
public class PropertyOption 
{
	@Id
	@GeneratedValue
	@Column(name = "poption_id")
	@XmlAttribute
	public Long id;

	@Column
	@XmlAttribute
	public String name;

	@ManyToOne
	@JoinColumn(name = "property_id")
	@XmlTransient
	public Property property;

	@XmlElement
	@Transient
	public Long countOfRecords; // For presentation purposes only

	@XmlElement
	@Transient
	public Long countExcluded;

	@XmlElement
	@Transient
	public Long countUniqueCompounds;

	@XmlAttribute
	@Transient
	public boolean multi; // means different options for one qualitative property/condition in batch editing

	public PropertyOption() {
	}

	public PropertyOption(String option, Property parent) {
		if(option == null || !validateOption(option))
			throw new UserFriendlyException("Property options should be a non-empty String and not a Float/Integer, but you provided: \"" + option + "\"");
		name = option.trim();
		property = parent;
	}

	public String toString()
	{
		return name;
	}

	@XmlAttribute
	public Boolean isReferenced()
	{
		if (ThreadScope.get().controller.equals("properties"))
		{
			Long count = (Long) Globals.session()
					.createQuery("select count(*) from ConditionSet cs join cs.values vs where vs.option=:option")
					.setParameter("option", this).uniqueResult();
			return count > 0;
		}
		return null;
	}

	public static PropertyOption getById(Long id)
	{
		return (PropertyOption) Globals.session().get(PropertyOption.class, id);
	}

	public boolean equals(Object obj)
	{
		return obj instanceof PropertyOption && ((PropertyOption) obj).name.equals(name);
	}

	public static boolean validateOption(String option) {
		String name = option.trim();
		if(name.length() == 0) return false;
		try {
			Double.valueOf(name);
			return false;
		}catch(NumberFormatException e){
		}
		return true;
	}
}
