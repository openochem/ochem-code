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

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Molecule;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.PropertyValue;
import qspr.util.MoleculePeer;
import qspr.util.UploadContext;
import qspr.util.WrapperThread;
import qspr.workflow.utils.QSPRConstants;

public class SolventConvertBasketRunner extends WrapperThread
{
	private static transient final Logger logger = LogManager.getLogger(SolventConvertBasketRunner.class);

	long basketId;

	public SolventConvertBasketRunner(long basket)
	{
		this.basketId = basket;
	}

	@SuppressWarnings("unchecked")
	@Override
	public void wrapped() throws Exception 
	{
		if(OCHEMConfiguration.autoLoginUser == null)return; // only for standalone

		Basket basket = Repository.basket.getById(basketId);
		Basket newBasket = Basket.getBasket(basket.name+"_solvent", Globals.userSession());

		UploadContext context = new UploadContext();

		List<Long> beIDs = Globals.session().createCriteria(BasketEntry.class).setProjection(Projections.groupProperty("id")).add(Restrictions.eq("basket", basket)).list();
		int batchSize = 50;
		int cnt = 0, dublicates = 0;

		Map<Integer,String> solvents = new HashMap<Integer,String>();

		Set<String> md5 = new HashSet<String>();
		for(BasketEntry entry: newBasket.entries)
			md5.add(entry.ep.md5); 

		Set<String> missed = new HashSet<String>();
		Set<String> missedO = new HashSet<String>();
		Set<ExperimentalProperty> missedP = new HashSet<ExperimentalProperty>();

		// Fetch and create the entries in batches
		for (int k = 0; k < beIDs.size(); k += batchSize)
		{
			if (cancelRequested)
				break;
			List<Long> batchIDs = beIDs.subList(k, Math.min(k + batchSize, beIDs.size()));
			List<BasketEntry> batchEntries = Globals.session().createCriteria(BasketEntry.class).add(Restrictions.in("id", batchIDs)).list();

			for (BasketEntry be : batchEntries)
			{
				if (cancelRequested)
					break;
				setStatus("Processed " + (cnt++) + " out of " + beIDs.size() + (dublicates > 0 ? "\n ("+dublicates+" dublicates skipped)" : ""));
				logger.info("Processed " + ((cnt)) + " out of " + beIDs.size() + (dublicates > 0 ? "\n ("+dublicates+" dublicates skipped)" : ""));
				ExperimentalProperty ep = be.ep;

				String mols[] = Various.molecule.splitOrderedByChargeAndSize(ep.molecule.getData());

				if(mols == null || mols.length == 0) {
					missedP.add(be.ep);
					continue;
				}

				ExperimentalProperty newEp = null;
				if(mols.length == 1) newEp = ep;
				else
				{

					newEp = be.ep.cloneForReference();
					newEp.article = ep.article;
					newEp.introducer = newEp.owner = Globals.userSession().user;
					newEp.rights = Globals.RIGHTS_NONE; // created by default as hidden
					newEp.connectedProperty = be.ep;
					newEp.molecule = MoleculePeer.fetchFromString(mols[0], context);

					newEp.moleculenames = null;

					String solvent = null;

					for(int i=1;i<mols.length;i++) {
						if(mols[i] == null) {
							missedP.add(be.ep);
							continue;
						}

						Molecule mol = MoleculePeer.fetchFromString(mols[i], context);
						String solv = solvents.containsKey(mol.mapping2.id)?solvents.get(mol.mapping2.id):Repository.record.getSolventName(mol.mapping2.id);
						if(solv == null) {
							if(Various.molecule.convertToCanonicalSMILES(mols[i]).equals("[I-]")) {
								mols[0] = Various.molecule.convertToCanonicalSMILES(mols[0])+"."+"[I-]";
								newEp.molecule = MoleculePeer.fetchFromString(mols[0], context);
							}else {
								missed.add(Various.molecule.convertToCanonicalSMILES(mols[i]));
								missedP.add(be.ep);
							}
							continue;
							//throw new UserFriendlyException("No solvent has been found for R" + ep.id + " SMILES: " + Various.molecule.convertToCanonicalSMILES(mols[i]));
						}
						solvents.put(mol.mapping2.id, solv);
						solvent = solvent == null? solv: solvent+";" + solv;
					}

					Property p = Property.getByName(QSPRConstants.SOLVENT_CONDITION);
					PropertyOption option = Repository.option.searchPropertyOptionByPermutations(solvent,p.id);
					if(option == null) {
						missedO.add(solvent);
						missedP.add(be.ep);
						continue;
						//throw new UserFriendlyException("No solvent option \"" + solvent +"\"was found for " + p.getName());
					}

					ConditionSet newCond =  new ConditionSet();
					newCond.values.add(new PropertyValue(option));
					if(newEp.conditions !=null)newEp.conditions.mergeAddOrUpdateWith(newCond);
					else
						newEp.conditions = newCond.get();

					newEp.updateHash();
					if (!newEp.hasConflicts())
						Globals.session().save(newEp);
					else
					{
						logger.info("A record already exists, using it.");
						dublicates++;
						if(md5.contains(newEp.md5))continue; // it was already stored 
						newEp = newEp.duplicate;
					}

				}

				if(!md5.contains(newEp.md5)) {
					BasketEntry newBe = new BasketEntry(newEp);
					newBe.basket = newBasket;
					newBe.ep = newEp;

					try{
						Globals.session().saveOrUpdate(newBe); // could be already there, just ignore exception
					}catch(Exception e){

					}
				}

				md5.add(newEp.md5);

			}

			//setStatus("Processed " + k + " out of " + basket.entries.size());
			Globals.restartAllTransactions(true);
		}

		newBasket = Basket.getBasket(basket.name+"missed_solvent", Globals.userSession());
		for(ExperimentalProperty ep: missedP) {
			BasketEntry newBe = new BasketEntry(ep);
			newBe.basket = newBasket;
			Globals.session().saveOrUpdate(newBe); // could be already there, just ignore exception
		}
		Globals.restartAllTransactions(true);

		System.out.println("solvents:");
		for(String s: missed) {
			System.out.println(s);
		}
		System.out.println("options:");
		for(String s: missedO
				) {
			System.out.println(s);
		}

		setStatus("Finished");
	}

}