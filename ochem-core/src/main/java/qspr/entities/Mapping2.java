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
import javax.persistence.Enumerated;
import javax.persistence.FetchType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.OneToMany;
import javax.persistence.OrderBy;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.workflow.utils.QSPRConstants;

@XmlRootElement
@Entity
public class Mapping2 
{
	@Id
	@Column(name = "mapping2_id")
	@javax.persistence.OrderBy
	@XmlAttribute
	public Integer id;

	@Column
	@XmlTransient
	public String inchi2;

	@ManyToOne //(fetch=FetchType.LAZY)
	@JoinColumn(name="mapping1_id", referencedColumnName="mapping1_id")
	public Mapping1 mapping1;

	@OneToMany(mappedBy = "mapping2", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@XmlTransient
	@OrderBy
	public List<Molecule> molecules = new ArrayList<Molecule>();

	@OneToMany(mappedBy = "mapping2", fetch = FetchType.LAZY)
	@XmlTransient
	public List<Mapping2Filter> filters = new ArrayList<Mapping2Filter>();

	@OneToMany(mappedBy = "mapping2", fetch = FetchType.LAZY)
	@XmlTransient
	public List<Mapping2Filter> basketFilters = new ArrayList<Mapping2Filter>();


	@OneToMany(mappedBy = "mapping", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	@XmlTransient
	@OrderBy
	public List<ValidatedFact> validatedFacts = new ArrayList<ValidatedFact>();

	/**
	 * 0 - not considered for MMP indexing,
	 * 1 - scheduled for indexing,
	 * 2 - submitted for indexing,
	 * 3 - indexed,
	 * 4 - error
	 */

	@Enumerated
	@Column
	public MMPAIndex mmpaIndexStatus = MMPAIndex.MMP_IGNORED; // by default and until manually requested

	public enum MMPAIndex{
		MMP_IGNORED, //(0)
		MMP_SCHEDULED,//(1)
		MMP_SUBMITTED,//(2)
		MMP_INDEXED,//(3)
		MMP_FAILED//(4)
	};
	/*
	public static final int MMP_IGNORED = 0;
	public static final int MMP_SCHEDULED = 1;
	public static final int MMP_SUBMITTED = 2;
	public static final int MMP_INDEXED = 3;
	public static final int MMP_FAILED = 4;
	 */

	@ManyToOne(fetch = FetchType.LAZY)
	@JoinColumn(name = "mapping2_id", insertable = false, updatable = false)
	@XmlTransient
	public StructureQuery xemistryIndex;

	public Molecule getMolecule() //Returning first (any) molecule for this mapping
	{
		Criteria criteria = Globals.session().createCriteria(Molecule.class)
				.add(Restrictions.eq("mapping2", this))
				.addOrder(Order.desc("id"))  // assuming that new data are better
				.setMaxResults(1);

		@SuppressWarnings("unchecked")
		List<Molecule> molecules = criteria.list();

		if (molecules.size() > 0)
			return molecules.get(0);
		return null;
	}

	public boolean isEmpty(){
		return inchi2 == null || inchi2.equals(QSPRConstants.EMPTY_MD5);
	}

	public boolean isBroken(){
		return inchi2 == null || inchi2.equals(QSPRConstants.EMPTY_MD5) || inchi2.length() != 25;
	}
}
