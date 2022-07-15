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

package qspr.modelling;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Mapping2;
import qspr.entities.Molecule;

public class ModelApplierAttachment implements Serializable
{
	private static final Logger logger = LogManager.getLogger(ModelApplierAttachment.class);
	private static final long serialVersionUID = 1L;
	
	public List<PredictedMolecule> molecules = new ArrayList<PredictedMolecule>();
	public Long basketId;
	public ModelApplierCache cache;
	
	public Basket getWorkData(boolean fetchMoleculesFromDb)
	{
		Basket basket = new Basket();
		basket.id = basketId;
		int i = 0;
		for (PredictedMolecule pMol : molecules)
		{
			ExperimentalProperty ep = new ExperimentalProperty();
			if (fetchMoleculesFromDb && pMol.epId != null)
			{
				ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, pMol.epId);
				if (i++ % 100 == 0)
					logger.info("Fetched " + i + " molecule structures from DB");
			}
			else
			{
				ep.id = pMol.epId;
				ep.molecule = new Molecule();
				ep.molecule.id = pMol.molId;
				ep.molecule.mapping2 = new Mapping2();
				ep.molecule.mapping2.id = pMol.mp2ID;
			}
			basket.addEntry(ep);
		}
		
		return basket;
	}
	
	public static void fetchMoleculesFromDatabase(Basket workdata)
	{
		logger.info("Fetching molecules from DB...");
		List<Long> idsToFetch = new ArrayList<Long>();

		for (int batchStart = 0, i = 0; i < workdata.entries.size(); i++)
		{
			BasketEntry be = workdata.entries.get(i);
			if (be.ep != null && be.ep.molecule != null && !be.ep.molecule.hasData())
				idsToFetch.add(be.ep.molecule.id);
			
			// Load records from DB in batches of 500 records. TODO: Consider restarting transaction (needs careful debugging) / Midnighter on Jun 27, 2011
			if (idsToFetch.size() >= 500 || i == workdata.entries.size() - 1)
			{
				if (!idsToFetch.isEmpty())
					Globals.session().createCriteria(Molecule.class).add(Restrictions.in("id", idsToFetch)).list();
				for (int k = batchStart; k <= i; k++)
				{
					BasketEntry be2 = workdata.entries.get(k);
					if (be2.ep != null && be2.ep.property == null && be2.ep.molecule.id > 0)
						be2.ep.molecule = Repository.molecule.getMolecule(be2.ep.molecule.id);
					batchStart = i + 1;
					idsToFetch.clear();
				}
				
				logger.info("" + i + " molecules fetched");
			}
		}
	}
	
	public ModelApplierAttachment setCache(ModelApplierCache cache)
	{
		this.cache = cache;
		return this;
	}
	
	public void setWorkData(Basket workdata)
	{
		basketId = workdata.id;
		for (BasketEntry be : workdata.entries)
			molecules.add(new PredictedMolecule(be));
	}
}