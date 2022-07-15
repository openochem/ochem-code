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

package qspr.fragmententities;

import java.io.Serializable;
import java.util.List;

import javax.persistence.Column;
import javax.persistence.Embeddable;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.IdClass;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Criteria;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;

@Entity
@XmlRootElement
@IdClass(Mapping1Fragment_Pk.class)
public class Mapping1Fragment {

	@XmlAttribute
	@Id
	Long mapping1_id;
	@XmlAttribute
	@Id
	Long fragment_id;

	@XmlTransient
	@ManyToOne
	@JoinColumn(name="mapping1_id", insertable=false, updatable=false)
	public Mapping1 mapping1;

	@XmlElement
	@ManyToOne
	@JoinColumn(name="fragment_id", insertable=false, updatable=false)
	public Fragment	fragment;

	@XmlElement
	@Column(name = "fragment_count")
	public Integer fragment_count;

	@SuppressWarnings("unchecked")
	public static Mapping1Fragment get(Mapping1 mapping, Fragment fragment, Integer count) {

		Criteria c = Globals.alternateSession().createCriteria(Mapping1Fragment.class)
				.add(Restrictions.eq("mapping1", mapping))
				.add(Restrictions.eq("fragment",fragment));

		List<Mapping1Fragment> list = c.list();
		Mapping1Fragment m1f;
		if (list.size() > 0)
			m1f = list.get(0);
		else
		{
			m1f = new Mapping1Fragment();
			m1f.fragment = fragment;
			m1f.mapping1 = mapping;
			m1f.fragment_id = fragment.id;
			m1f.mapping1_id = mapping.id;
		}

		m1f.fragment_count = count;

		Globals.alternateSession().saveOrUpdate(m1f);
		return m1f;
	}

	@SuppressWarnings("unchecked")
	public static List<Long> getMapping1IdsForFragment(Long fragment_id) {

		Criteria c = Globals.alternateSession().createCriteria(Mapping1Fragment.class)
				.setProjection(Projections.property("mapping1_id"))
				.add(Restrictions.eq("fragment_id", fragment_id))
				.setMaxResults(10000);

		List<Long> list = c.list();


		return list;
	}
}

@Embeddable
class Mapping1Fragment_Pk implements Serializable
{
	private static final long serialVersionUID = 1L;
	Long mapping1_id;
	Long fragment_id;

	public boolean equals(Mapping1Fragment_Pk mp1f_pk)
	{
		if (mapping1_id == null || fragment_id == null)
			return false;

		if (mapping1_id.equals(mp1f_pk.mapping1_id) && fragment_id.equals(mp1f_pk.fragment_id))
			return true;

		return false;
	}

}