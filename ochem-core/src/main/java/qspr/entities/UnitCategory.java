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

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.OneToMany;
import javax.persistence.OrderBy;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.ThreadScope;

/**
 * Category of measurement units (e.g., concentration, time, mass)
 */
@Entity
@XmlRootElement(name = "unitcategory")
public class UnitCategory {

	@Id
	@GeneratedValue
	@Column(name = "category_id")
	@XmlAttribute
	public Long id;

	@Column
	@XmlAttribute
	public String name;

	@OneToMany(mappedBy = "unitCategory")
	@XmlTransient
	public Set<Property> properties = new HashSet<Property>();

	@OneToMany(mappedBy = "category")
	//@XmlElement(name = "unit")
	@XmlTransient
	@OrderBy("name")
	public Set<Unit> units;// = new HashSet<Unit>();

	@Transient
	@XmlAttribute
	public Boolean selected;

	@XmlElement(name = "unit")
	private Set<Unit> getUnits()
	{
		if (Globals.getMarshallingOption(MarshallingOption.UNITCATEGORY_UNITS))
			return units;
		else
			return null;
	}

	public String toString()
	{
		return this.name;
	}

	@XmlTransient@Transient
	private Unit defaultUnit;

	@XmlTransient@Transient
	public Unit getDefaultUnit()
	{
		if (defaultUnit != null)
			return defaultUnit;

		Criteria c = Globals.session().createCriteria(Unit.class)
				.add(Restrictions.eq("isdefault", 1))
				.add(Restrictions.eq("category", this));
		@SuppressWarnings("unchecked")
		List<Unit> list = c.list();

		if (list.size() > 0)
			return defaultUnit = list.get(0);

		if (units != null && units.size() > 0)
		{
			Object[] units = this.units.toArray();
			return defaultUnit = (Unit)units[0];
		}

		return null;
	}

	@Transient@XmlAttribute(name = "default-unit")
	public String getDefaultUnitName()
	{
		if ("unit".equals(ThreadScope.get().controller) || "properties".equals(ThreadScope.get().controller))
			return getDefaultUnit().getName();
		return null;
	}

	//	@XmlElement(name = "unit")
	//	public List<Unit> getUnits()
	//	{
	//		if (null != units)
	//			Collections.sort(units);
	//		return units;
	//	}

	public boolean equals(Object cat)
	{
		if (cat != null && cat instanceof UnitCategory)
			return id.equals(((UnitCategory)cat).id);
		return false;
	}

	public Unit getUnitByName(String name)
	{
		return (Unit) Globals.session().createCriteria(Unit.class).add(Restrictions.eq("name", name)).add(Restrictions.eq("category", this)).uniqueResult();
	}

	public static UnitCategory getByName(String name)
	{
		@SuppressWarnings("unchecked")
		List<UnitCategory> l = Globals.session().createCriteria(UnitCategory.class).add(Restrictions.like("name",name)).list();
		if (l.size() > 0)
			return l.get(0);
		UnitCategory uc = new UnitCategory();
		uc.name = name;
		return uc;
	}
}
