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

import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.OrderBy;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.workflow.utils.QSPRConstants;


@Entity
@XmlRootElement(name = "ValidatedFact")
public class ValidatedFact // the class name is misleading. Something like "ValidatedName" would be better
{

	public static final int SOURCE_TOXICITY = 0;
	public static final int SOURCE_PUBCHEM = 1;
	public static final int VALIDATED = 1;
	public static final int SYNONYM = 0;

	@Id
	@Column(name="validatedfact_id")
	@GeneratedValue
	@OrderBy
	@XmlAttribute
	public Long id;

	@ManyToOne
	@JoinColumn(name = "mapping2_id")
	public Mapping2 mapping;

	@ManyToOne
	@JoinColumn(name = "molecule_name_id", referencedColumnName="molecule_name_id")
	public MoleculeName moleculename;

	@Column(name = "source")
	public int source;

	@Column(name = "source_id")
	public int sourceid;

	@Column(name = QSPRConstants.VALIDATION)
	public int validated;

	public ValidatedFact()
	{

	}

	@SuppressWarnings("unchecked")
	public static ValidatedFact getPrimary(MoleculeName name)
	{
		List<ValidatedFact> vFacts = Globals.session().createCriteria(ValidatedFact.class)
				.add(Restrictions.eq("moleculename", name))
				.add(Restrictions.eq("validated", ValidatedFact.VALIDATED))
				.setMaxResults(1)
				.list();

		if (vFacts.size() > 0)
			return vFacts.get(0);
		else
			return null;
	}

	public static ValidatedFact getSynonym(MoleculeName name, Mapping2 mol)
	{
		Criteria c = Globals.session().createCriteria(ValidatedFact.class)
				.add(Restrictions.eq("moleculename", name))
				.add(Restrictions.eq("validated", ValidatedFact.SYNONYM))
				.setMaxResults(1);

		if (mol == null)
			c.add(Restrictions.isNull("mapping"));
		else
			c.add(Restrictions.eq("mapping", mol));

		@SuppressWarnings("unchecked")
		List<ValidatedFact> vFacts = c.list();

		if (vFacts.size() > 0)
			return vFacts.get(0);
		else
			return null;
	}
	//	public static ValidatedFact get(MoleculeName name, Mapping2 mapping)
	//	{
	//		Criteria c = Globals.session().createCriteria(ValidatedFact.class);
	//			c.add(Restrictions.eq("moleculename", name));
	//			if (mapping == null)
	//				c.add(Restrictions.isNull("mapping"));
	//			else
	//				c.add(Restrictions.eq("mapping", mapping));	 
	//			
	//			List list = c.list();
	//		
	//		if (list.size() > 0)
	//			return (ValidatedFact) list.get(0);
	//		else
	//		{
	//			return null;
	//		}
	//	}
}
