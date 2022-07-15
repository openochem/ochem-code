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

package com.eadmet.business;

import java.util.ArrayList;
import java.util.List;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.AbstractMoleculeName;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Molecule;
import qspr.entities.MoleculeName;
import qspr.entities.ValidatedFact;
import qspr.util.NCBI_Utility;

public class MoleculeNamesService 
{
	public ExperimentalProperty checkNames(ExperimentalProperty ep) 
	{
		for (MoleculeName molName : ep.moleculenames) 
		{
			boolean validated = false;
			for (ValidatedFact vf : molName.validatedFacts) 
			{
				if (vf.validated == ValidatedFact.VALIDATED) 
				{
					validated = true;
					break;
				}
			}

			if (!validated) 
			{
				try 
				{
					NCBI_Utility.getMoleculeByName(molName.name.trim());
				} catch (Exception e) 
				{
					e.printStackTrace();
				}
			}
		}
		return ep;
	}

	public ValidatedFact validate(MoleculeName mn, Molecule m) 
	{
		ValidatedFact vf = ValidatedFact.getPrimary(mn);

		if (vf == null)
			vf = new ValidatedFact();

		vf.mapping = m.mapping2;
		vf.moleculename = mn;
		vf.source = ValidatedFact.SOURCE_TOXICITY;
		vf.sourceid = Globals.userSession().user.id.intValue();
		vf.validated = ValidatedFact.VALIDATED;
		Globals.session().saveOrUpdate(vf);
		return vf;
	}

	public void invalidate(MoleculeName mn) 
	{
		@SuppressWarnings("unchecked")
		List<ValidatedFact> vFacts = Globals.session().createCriteria(ValidatedFact.class)
				.add(Restrictions.eq("moleculename", mn))
				.add(Restrictions.eq("validated", ValidatedFact.VALIDATED))
				.setMaxResults(1)
				.list();

		ValidatedFact vf;
		if (vFacts.size() > 0) 
		{
			vf = vFacts.get(0);
			Globals.session().delete(vf);
			if (mn.validatedFacts.contains(vf))
				mn.validatedFacts.remove(vf);
		}

	}

	public void setNamesForEP(ExperimentalProperty ep, String[] names)
	{
		List<MoleculeName> nameList = new ArrayList<MoleculeName>();
		if (names != null)
		{
			for (int i = 0; i < names.length; i++)
			{
				String name = names[i].trim();
				MoleculeName currentMolName = ep.containsName(name); 
				if (currentMolName == null)
				{
					MoleculeName molname = MoleculeName.get(name);
					if (molname != null)
						nameList.add(molname);
				} else 
				{
					nameList.add(currentMolName);
				}
			}
		}

		ep.moleculenames.clear();
		ep.moleculenames.addAll(nameList);
		Globals.session().saveOrUpdate(ep);
	}

	public List<AbstractMoleculeName> listMoleculeNames(ExperimentalProperty ep, Molecule mol)
	{
		List<AbstractMoleculeName> moleculenames;
		if (mol != null) 
		{
			ep.molecule = mol;
			moleculenames = ep.getCheckedNamesIndividual(); //Here we get names "old way" - by coloring each individually
		} else
			moleculenames = ep.getCheckedNames();		
		return moleculenames;
	}

	public List<MoleculeName> listSynonyms(Molecule m, ExperimentalProperty ep, PaginationFilter pager)
	{
		if(!m.mapping2.isEmpty())
		{
			ResultBuilder<MoleculeName> rb = getResultBuilder(m, ep).withPager(pager);
			pager.totalSize = rb.count();
			return rb.list();
		} else
		{
			pager.totalSize = 0L;
			return new ArrayList<MoleculeName>();
		}
	}

	public List<MoleculeName> listSynonymsFromSearch(Molecule m)
	{
		List<MoleculeName> names = new ArrayList<MoleculeName>();
		if (m != null && m.searchSynonymList != null)
		{
			for (String synonym : m.searchSynonymList) 
			{
				MoleculeName molName = MoleculeName.get(synonym);
				if (molName != null)
					names.add(molName);
			}
		}
		return names;
	}

	private ResultBuilder<MoleculeName> getResultBuilder(Molecule m, ExperimentalProperty ep)
	{
		Criteria c = Globals.session()
				.createCriteria(Molecule.class)
				.add(Restrictions.eq("id", m.id))
				.createCriteria("mapping2")
				.createCriteria("molecules")
				.createCriteria("experimentalProperties");

		if (ep != null)
			c = c.add(Restrictions.ne("id", ep.id));

		c = c.createCriteria("moleculenames");

		return new ResultBuilder<MoleculeName>(c, MoleculeName.class, Globals.session()).setDistinct(true); 
	}
}
