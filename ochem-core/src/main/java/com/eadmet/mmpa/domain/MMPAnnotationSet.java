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
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Article;
import qspr.entities.Property;
import qspr.entities.User;

@Entity
@XmlRootElement
public class MMPAnnotationSet
{

	@Id
	@GeneratedValue
	@Column(name = "as_id")
	@XmlAttribute
	public Long id;

	@ManyToOne
	@JoinColumn(name = "user_id")
	public User user;

	@Column
	public String name;

	// For marshalling to XML
	@Transient
	public List<Property> properties;

	@Column
	public boolean published;

	@Column
	public boolean featured;

	@Column
	public String description;

	@XmlTransient
	@ManyToOne
	@JoinColumn(name = "article_id")
	public Article article;

	public List<Long> getIdList()
	{
		List<Long> l = new ArrayList<Long>();
		l.add(id);
		return l;
	}

	@SuppressWarnings("unchecked")
	public static MMPAnnotationSet getByName(String name)
	{
		Criteria c = Globals.session().createCriteria(MMPAnnotationSet.class).add(Restrictions.like("name", name));
		if (Globals.userSession().user != null)
			c.add(Restrictions.eq("user", Globals.userSession().user));
		List<MMPAnnotationSet> list = c.list();
		if (list.size() > 0)
			return list.get(0);
		MMPAnnotationSet set = new MMPAnnotationSet();
		set.name = name;
		set.user = Globals.userSession().user;
		Globals.session().saveOrUpdate(set);
		return set;
	}

	public static MMPAnnotationSet getById(Long id) {
		return (MMPAnnotationSet) Globals.session().get(MMPAnnotationSet.class, id);
	}
}
