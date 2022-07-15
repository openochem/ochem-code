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

import java.io.IOException;
import java.util.List;
import java.util.Set;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Molecule;

@XmlRootElement(name = "fragment")
public class Fragment {

	@XmlAttribute
	public Long id;

	@XmlElement
	public String fragment_data;

	@XmlElement
	public String inchi1;

	@XmlTransient
	public Set<Mapping1Fragment> mapping1Fragments;

	@XmlElement
	public Long size;

	public String getFragment() {
		return fragment_data;
	}

	@Override // needed for equality in a set
	public boolean equals(Object obj)
	{
		Fragment other = (Fragment) obj;
		return inchi1.equals(other.inchi1);
	}

	@Override // needed for equality in a set
	public int hashCode()
	{
		return (this.inchi1).hashCode();
	}

	public static Fragment get(String inchi1, String data) {

		Fragment f = get(inchi1); 
		// if this inchi is not there create a new fragment, add it to db and return it
		if (f == null) 
		{
			Fragment newFragment = new Fragment();
			newFragment.fragment_data = data;
			newFragment.inchi1 = inchi1;
			Globals.alternateSession().save(newFragment);
			return newFragment;
		} else 
		{
			return f;
		}
	}

	@SuppressWarnings("unchecked")
	public static Fragment get(String inchi1) {
		Criteria criteria = Globals.alternateSession().createCriteria(Fragment.class)
				.add(Restrictions.eq("inchi1", inchi1));

		List<Fragment> fragment = criteria.list();

		if (fragment.size() > 0) 
			return fragment.get(0);
		else 
			return null; 

	}

	public String[] calInchi(String _data) throws IOException, InterruptedException{
		return Molecule.getInChiKeys(_data);
	}
}
