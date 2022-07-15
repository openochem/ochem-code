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

import javax.persistence.Embeddable;
import javax.persistence.Entity;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;

@Entity
@XmlRootElement
public class FragmentTree implements Serializable{

	//	/**
	//	 * 
	//	 */
	private static final long serialVersionUID = 1342966780843953968L;
	//
	@XmlAttribute
	@Id
	Long fragment_id;

	@XmlAttribute
	@Id
	Long child1_id;

	@XmlAttribute
	Long child2_id;

	@XmlElement
	@ManyToOne
	@JoinColumn(name="fragment_id", insertable=false, updatable=false)
	public Fragment	fragment;

	@XmlElement
	@ManyToOne
	@JoinColumn(name="child1_id", insertable=false, updatable=false)
	public Fragment	child1;

	@XmlElement
	@ManyToOne
	@JoinColumn(name="child2_id", insertable=false, updatable=false)
	public Fragment	child2;

	@SuppressWarnings("unchecked")
	public static FragmentTree get(Fragment fragment, Fragment child1, Fragment child2) {

		Criteria c = Globals.alternateSession().createCriteria(FragmentTree.class)
				.add(Restrictions.eq("fragment", fragment))
				.add(Restrictions.eq("child1", child1));

		List<FragmentTree> list = c.list();
		FragmentTree ft;
		if (list.size() > 0)
			ft = list.get(0);
		else
		{
			ft = new FragmentTree();
			ft.fragment = fragment;
			ft.fragment_id = fragment.id;

			ft.child1 = child1;
			ft.child2 = child2;

			//			ft.child1_id = (ft.child1 != null) ? child1.id : Long.valueOf(0);
			//			ft.child2_id = (ft.child1 != null) ? child2.id : Long.valueOf(0);

			ft.child1_id = child1.id;
			ft.child2_id = child2.id;

		}

		Globals.alternateSession().saveOrUpdate(ft);
		return ft;
	}

}

@Embeddable
class FragmentTree_Pk implements Serializable
{
	private static final long serialVersionUID = 1L;
	Long fragment_id;
	Long child1_id;

	public boolean equals(FragmentTree_Pk ft_pk)
	{
		if (fragment_id == null || child1_id == null)
			return false;

		if (fragment_id.equals(ft_pk.fragment_id) && child1_id.equals(ft_pk.child1_id))
			return true;

		return false;
	}

}