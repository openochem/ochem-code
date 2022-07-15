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

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.JoinTable;
import javax.persistence.ManyToMany;
import javax.persistence.OneToMany;
import javax.persistence.OrderBy;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;

import com.eadmet.utils.OCHEMUtils;

@Entity
@XmlRootElement(name = "moleculename")
public class MoleculeName extends AbstractMoleculeName
{
	@Id
	@GeneratedValue
	@Column(name = "molecule_name_id")
	@XmlAttribute
	public Long id;

	@XmlAttribute
	public String name;

	@XmlAttribute
	public String md5;

	@Transient
	@XmlAttribute
	public int validation;

	@Transient
	@XmlAttribute
	public String user;


	@OneToMany(mappedBy = "moleculename", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@XmlTransient
	@OrderBy
	public List<ValidatedFact> validatedFacts = new ArrayList<ValidatedFact>();


	@ManyToMany
	(
			targetEntity = ExperimentalProperty.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE}
			)
	@JoinTable
	(
			name="ExperimentalPropertyName",
			inverseJoinColumns={@JoinColumn(name="exp_property_id")},
			joinColumns={@JoinColumn(name="molecule_name_id")}
			)
	@XmlTransient
	public List<ExperimentalProperty> records = new ArrayList<ExperimentalProperty>();

	//	@ManyToMany(mappedBy = "moleculenames")
	//	Set<MoleculeName> names = new HashSet<MoleculeName>();

	public static String nameMD5(String name)
	{
		return OCHEMUtils.getMD5(unifyName(name));
	}

	public static String unifyName(String name)
	{
		return name.toLowerCase();
	}

	public static MoleculeName get(String name)
	{
		if (name == null || "".equals(name))
			return null;

		@SuppressWarnings("rawtypes")
		List mNames = Globals.session().createCriteria(MoleculeName.class).add(Restrictions.eq("md5", nameMD5(name))).list();

		if (mNames.size() > 0)
			return (MoleculeName) mNames.get(0);
		else
		{
			MoleculeName mName = new MoleculeName();
			mName.name = name;
			mName.md5 = nameMD5(name);
			Globals.session().saveOrUpdate(mName);
			return mName;
		}
	}

	@Override
	public boolean equals(Object o)
	{
		MoleculeName mn = (MoleculeName)o;

		if (id == null)
			return false;

		return (id.equals(mn.id));
	}


	@Override
	public int hashCode()
	{
		return id.hashCode();
	}

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
		return validation;
	}
}
