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

import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.JoinTable;
import javax.persistence.ManyToMany;
import javax.persistence.OneToMany;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Criteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.dao.Various;

@Entity
@XmlRootElement(name = "mapping")
public class Mapping1 
{
	@Id
	@Column(name = "mapping1_id")
	@javax.persistence.OrderBy
	@XmlAttribute
	public Long id;

	@Column
	@XmlTransient
	public String inchi1;

	@Transient
	public Boolean isFragmented;

	// Some confidential molecules that are not referenced in any tags or records should not be shown to the user
	// This is controlled by this flag, which should be regularly (TODO) updated in some cron job / Midnighter on Jun 28, 2011
	@Column@XmlTransient
	public boolean visible;

	@XmlAttribute
	public boolean isSelected()
	{
		if (Globals.userSession() == null)
			return false;
		Set<Long> selectionList = Globals.userSession().selectionMoleculeList;
		return selectionList.contains(this.id);
	}	
	@OneToMany(mappedBy = "mapping1", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@XmlTransient
	public List<qspr.entities.Molecule> molecules = new ArrayList<qspr.entities.Molecule>();

	@OneToMany(mappedBy = "mapping1", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@XmlTransient
	public List<qspr.entities.Mapping2> mapping2 = new ArrayList<qspr.entities.Mapping2>();


	@ManyToMany
	(
			targetEntity = Fragment.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE}
			)
	@JoinTable
	(
			name = "Mapping1Fragment",
			joinColumns = {@JoinColumn(name="mapping1_id")},
			inverseJoinColumns ={@JoinColumn(name="fragment_id")}
			)
	@XmlTransient
	public Set<Fragment> moleculeFragment = new HashSet<Fragment>();


	@Transient
	qspr.entities.Molecule molecule;

	@ManyToMany
	(
			targetEntity = Tag.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE},
			fetch = FetchType.LAZY
			)
	@JoinTable
	(
			name="MoleculeTag",
			joinColumns={@JoinColumn(name="mapping1_id")},
			inverseJoinColumns={@JoinColumn(name="tag_id")}
			)
	@XmlTransient
	public Set<Tag> tags = new HashSet<Tag>();	

	@SuppressWarnings("unchecked")
	@XmlElement(name = "tag")
	public List<Tag> getTags()
	{
		if ("molbrowser".equals(ThreadScope.get().controller))
		{
			Disjunction rights = Restrictions.disjunction();
			rights.add(Restrictions.eq("isPublic", Boolean.TRUE));
			if (Globals.userSession().user != null)
				rights.add(Restrictions.eq("introducer", Globals.userSession().user));

			return Globals.session().createCriteria(Tag.class)
					.createAlias("mapping", "mp")
					.add(Restrictions.eq("mp.id", this.id))
					.add(rights).list();
		}

		return null;
	}


	public static Mapping1 get(String inchi1, String data)
	{
		Criteria criteria = Globals.session().createCriteria(Mapping1.class)
				.add(Restrictions.eq("inchi1", inchi1));

		@SuppressWarnings("unchecked")
		List<Mapping1> mapping1s = criteria.list();

		if (mapping1s.size() == 0)
		{
			Mapping1 mapping1 = new Mapping1();
			mapping1.inchi1 = inchi1;

			qspr.fragmententities.Mapping1 externalMapping = qspr.fragmententities.Mapping1.get(inchi1, data);			//Get Mapping1 ID from Central Repository...
			mapping1.id = externalMapping.id;     

			Globals.session().saveOrUpdate(mapping1);
			return mapping1;
		}
		else
			return mapping1s.get(0);
	}

	@XmlElement(name = "exp-property")
	@SuppressWarnings("unchecked")
	public List<Object> getExperimentalProperty()
	{
		if (ThreadScope.get().controller.equals("molbrowser"))
		{
			Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class)
					.add(Restrictions.isNull("deleted"));
			ExperimentalProperty.addAccessRestrictions(criteria, Globals.RIGHTS_FREELY_AVAILABLE, Globals.userSession().user, false, true);

			criteria.createAlias("molecule", "molecule")
			.add(Restrictions.eq("molecule.mapping1", this)).setMaxResults(6).addOrder(Order.asc("property"));
			return criteria.list();
		}
		else
			return null;
	}

	@XmlAttribute(name="mol_id")
	public Long getMolId()
	{
		if (ThreadScope.get().controller.equals("molbrowser"))
		{
			getMolecule();
			if(molecule != null)
				return molecule.id;
		}
		return null;
	}

	@XmlAttribute(name="mp2")
	public Integer getMp2()
	{
		if (ThreadScope.get().controller.equals("molbrowser"))
		{
			getMolecule();
			if(molecule != null)
				return molecule.mapping2.id;
		}
		return null;
	}

	@XmlAttribute(name="MF")
	public String getMF() throws IOException
	{
		if (ThreadScope.get().controller.equals("molbrowser"))
		{
			return Various.molecule.getFormula(getMolecule());
		}
		return null;
	}

	@XmlElement(name="MW")
	public String getMW() throws IOException
	{
		if (ThreadScope.get().controller.equals("molbrowser"))
		{
			DecimalFormat format = new DecimalFormat ( ".###" ) ;
			return format.format(Various.molecule.getMass(getMolecule()));
		}
		return null;
	}

	@XmlElement(name="smiles")
	public String getSmiles()
	{
		if (ThreadScope.get().controller.equals("molbrowser") && Globals.isValidatedUser())
			try {

				return Various.molecule.convertToSmilesOrSmart(getMolecule());
			} catch (IOException e) {
				e.printStackTrace();
				return  "no smiles";
			} 
		return null;
	}

	@XmlElement(name="inchikey")
	public String getInchiKey()
	{
		if (ThreadScope.get().controller.equals("molbrowser"))
		{
			getMolecule();
			if(this.inchi1.length() == 14)
				return molecule.mapping2.inchi2;
			else
				return "InChIKey can not be calculated";
		}
		return null;
	}


	@XmlTransient
	private String getMolecule()
	{
		if (molecule == null)
		{
			Criteria criteria = 
					Globals.session().createCriteria(qspr.entities.Molecule.class)
					.add(Restrictions.eq("mapping1", this));
			if(criteria.list().size() > 0)
				molecule = (qspr.entities.Molecule)criteria.list().get(0);
		}

		return molecule.getData();
	}

	public static void updateVisibilityFlag()
	{
		// TODO: This should go in some cron job
		Globals.session().createSQLQuery("update Mapping1, MoleculeTag set visible=1 where (visible=0) and (Mapping1.mapping1_id=MoleculeTag.mapping1_id)").executeUpdate();
		Globals.restartMainTransaction(true);
	}
}
