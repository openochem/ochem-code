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

import java.util.ArrayList;
import java.util.List;

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

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.ProjectionList;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.ThreadScope;
import qspr.annotations.Loggable;
import qspr.util.unitconversion.UnitConversion;

import com.eadmet.exceptions.UserFriendlyException;

/**
 * Unit of measurement (e.g., kg, centimeter, second)
 */
@Entity
//@Loggable
@XmlRootElement(name = "unit")
public class Unit implements Comparable<Unit>
{
	private static transient final Logger logger = LogManager.getLogger(Unit.class);

	@Id
	@GeneratedValue
	@Column(name = "unit_id")
	@XmlAttribute
	public Long id;

	private String name;

	/**
	 * Used by XML
	 */

	public String shortName;

	@Column
	@XmlElement
	public String description;

	@Column
	@XmlAttribute
	public int isdefault;

	@Column(name = "to_default_conversion")
	@Loggable
	public String toDefaultConversion;

	@Column(name = "from_default_conversion")
	@Loggable
	public String fromDefaultConversion;

	@ManyToOne
	@JoinColumn(name = "category_id", referencedColumnName="category_id")
	@XmlTransient
	@Loggable
	public UnitCategory category;

	@ManyToOne
	@JoinColumn(name = "modifier_id")
	@XmlTransient
	@Loggable(name = "modifier")
	public User owner;

	@ManyToOne
	@JoinColumn(name = "introducer_id")
	@XmlTransient
	@Loggable
	public User introducer;

	@Column(name = "doc_term")
	@XmlElement
	public String documentationTerm;

	@XmlElement
	@Transient
	public Long countInProperty; // Solely for marshalling purpose.. 

	@XmlAttribute
	@Transient
	public boolean multi; // means different units for one quantitative property / condition in batch editing

	public String toString()
	{
		return this.name;
	}

	public void doCheckRights() throws Exception
	{
		if (!(this.owner == null || this.owner.equals(Globals.userSession().user)))
			if (!(Globals.userSession().user != null && Globals.userSession().user.rank > this.owner.rank))
			{
				if (!(Globals.userSession().user != null && Globals.userSession().user.rank.equals(this.owner.rank)))
				{
					logger.info("Not permitted, delete has been ignored");
					throw new UserFriendlyException("Not permitted, delete has been ignored");
				}
			}
	}

	@XmlElement(name = "unitCategory")
	private UnitCategory getUnitCategory()
	{
		if (ThreadScope.get().controller.equals("unit"))
			return category;
		else
			return null;
	}

	@XmlAttribute(name="unit-record")
	public Long getPropertyCount()
	{
		if (Globals.getMarshallingOption(MarshallingOption.UNIT_RECORDS_COUNT))
		{
			Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class)
					.add(Restrictions.isNull("deleted"))
					.add(Restrictions.eq("unit", this));
			//.createAlias("conditions", "c").createAlias("c.values", "cv")
			//.add(Restrictions.or(Restrictions.eq("unit", this), Restrictions.eq("cv.unit", this)));
			ExperimentalProperty.addAccessRestrictions(criteria, Globals.RIGHTS_FREELY_AVAILABLE, Globals.userSession().user, false, true);

			criteria.setProjection(Projections.countDistinct("id"));
			return (Long)criteria.list().get(0);
		}
		return 0L;
	}

	@XmlElement(name = "unit-property")
	private List<Property> getProperties()
	{
		if (Globals.getMarshallingOption(MarshallingOption.UNIT_PROPERTIES))
		{
			if (id == null)
				return null;

			List<Property> properties = new ArrayList<Property>();

			ProjectionList projList = Projections.projectionList();
			projList.add(Projections.groupProperty("property"));
			projList.add(Projections.countDistinct("id"));

			@SuppressWarnings("unchecked")
			List<Object[]> res = Globals.session().createCriteria(ExperimentalProperty.class)
			.add(Restrictions.eq("unit", this))
			.setProjection(projList).list();

			for (Object[] objects : res) 
			{
				Property property = (Property) objects[0];
				Property copy = new Property();
				copy.id = property.id;
				copy.setName(property.getName());
				copy.count = Long.valueOf((Long)objects[1]);
				properties.add(copy);
			}

			return properties;
		}
		return null;
	}

	@XmlElement(name = "owner")
	public String getOwner()
	{
		if (owner != null)
			return owner.login;
		else
			return null;
	}

	@XmlElement(name = "introducer")
	public String getIntroducer()
	{
		if (introducer != null)
			return introducer.login;
		else
			return null;
	}

	@XmlAttribute
	public Unit setName(String _name)
	{
		name = _name.trim();
		shortName = Unit.shortName(name);
		return this;
	}

	public String getName()
	{
		return name;
	}

	public static String shortName(String name)
	{
		String result = name.toLowerCase();//.replaceAll("[^0-9a-z\\+\\-\\%]*", "");
		return result;
	}

	public int compareTo(Unit o)
	{
		return name.toLowerCase().compareTo(o.getName().toLowerCase());
	}

	public boolean equals(Unit otherUnit)
	{
		return otherUnit != null && otherUnit.id != null && otherUnit.id.equals(id);
	}

	public boolean isOppositeTo(Unit anotherUnit)
	{
		if (this.equals(anotherUnit))
			return false;
		// Check if increasing the value in this unit leads to decreasing in the other unit
		double val1 = UnitConversion.convert(1.0, this, anotherUnit, 100.0);
		double val2 = UnitConversion.convert(1.1, this, anotherUnit, 100.0);
		return val2 < val1;
	}

	public static Unit getByNameAndCategory(String name, String category, boolean createIfMissing)
	{
		@SuppressWarnings("unchecked")
		List<Unit> l =	Globals.session().createCriteria(Unit.class)
		.add(Restrictions.eq("shortName", Unit.shortName(name)))
		.createAlias("category", "c")
		.add(Restrictions.eq("c.name", category))
		.list();

		if (l.size() > 0)
			return l.get(0);

		if (!createIfMissing)
			return null;

		Unit u = new Unit();
		u.setName(name);
		u.category = UnitCategory.getByName(category);
		return u;
	}

	public static Unit getById(Long id)
	{
		return (Unit) Globals.session().get(Unit.class, id);
	}
}
