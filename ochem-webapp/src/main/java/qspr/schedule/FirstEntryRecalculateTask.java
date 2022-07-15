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

package qspr.schedule;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.dao.Repository;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Property;

// Experimental feature "Primary records"

@DatabaseMaintenanceJob
public class FirstEntryRecalculateTask extends OchemCronjobTask {

	public static void main(String[] args) 
	{
		new FirstEntryRecalculateTask().executeTask();
	}

	@SuppressWarnings("unchecked")
	private void resetReferencingRecords(ExperimentalProperty ep)
	{
		try 
		{
			List<ExperimentalProperty> referencingRecords = 
					Globals.session().createCriteria(ExperimentalProperty.class).add(Restrictions.or(
							Restrictions.eq("firstEntry", ep.id))).list();

			if(referencingRecords.size() > 0)
				log(referencingRecords.size() + " referencing records found for record R" + ep.id);

			for (ExperimentalProperty referencingRecord : referencingRecords) 
			{
				referencingRecord.firstEntry = null;
				log("Invalidating a referenced record " + referencingRecord.id);
			}

		} catch (Exception e) {
			log("ERROR while fetching referencing records for ep id " + ep.id);
			e.printStackTrace();
		}
	}

	@SuppressWarnings("unchecked")
	public void executeTask()
	{	

		ThreadScope.get().disableTrackChanges = true;
		long timer = System.nanoTime();
		Globals.startMainTransaction();

		List<Long> ids = new ArrayList<Long>();


		//ids = Repository.basket.getById(1373l).getIds();
		//ids = new ArrayList<Long>(); ids.add(51207840l);

		if(ids.size() == 0)
			ids = Globals.session().createCriteria(ExperimentalProperty.class)
			.add(Restrictions.isNull("firstEntry"))
			.add(Restrictions.isNull("deleted"))
			//.createAlias("owner", "o")
			//.createAlias("introducer", "i")
			//.createAlias("article", "a")
			//.add(Restrictions.not(Restrictions.like("o.login", OCHEMTestHelpers.TEST_USER_PREFIX+"%")))
			//.add(Restrictions.not(Restrictions.like("i.login", OCHEMTestHelpers.TEST_USER_PREFIX+"%")))
			.setProjection(Projections.id())
			.list();

		Set<Long> recordIdSet = new HashSet<Long>();
		Collections.sort(ids);
		recordIdSet.addAll(ids);

		int lastRestartSize, totalSize;
		lastRestartSize = totalSize = ids.size();

		logger.info("Got Total number of records: " + totalSize);

		while (ids.size() > 0)
		{
			Long recordId = ids.remove(0); // at least one is removed
			ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, recordId);

			if(ep == null) continue; // should not happen unless deleted during same update

			if (ep.article.publicationDate == null)
				ep.article.publicationDate = new Date();

			if(ep.isDummyOrEmpty()) { // we just skip such records
				if(ep.firstEntry == null || (long)ep.firstEntry != (long)ep.id) {
					ep.firstEntry = ep.id;
					Globals.session().saveOrUpdate(ep);
				}
				continue;
			}

			List<ExperimentalProperty> similarRecords = null;


			Property p = Repository.property.getPropertyById(ep.property.id);
			Map<String,List<ExperimentalProperty>> map = new HashMap<String,List<ExperimentalProperty>>();

			resetReferencingRecords(ep);
			similarRecords = getSimilarRecords(ep);

			for(ExperimentalProperty eps:similarRecords) {
				String hash = eps.getObligatoryConditionHash(p); 
				if(!map.containsKey(hash))
					map.put(hash, new ArrayList<ExperimentalProperty>());
				List<ExperimentalProperty> obl = map.get(hash);
				obl.add(eps);
			}

			if(similarRecords == null || similarRecords.size() < 2){ // majority of records will have only itself
				ep.firstEntry = ep.id;
				Globals.session().saveOrUpdate(ep);
			}else {
				logger.info("Found for record: " + ep.id + " in total  " + similarRecords.size() + " separated in " + map.size());

				for(String hash: map.keySet()){
					similarRecords = map.get(hash);

					ExperimentalProperty primaryRecord = null;

					// we look for the primary record for all records, including our record
					for (ExperimentalProperty record : similarRecords)
					{
						if(primaryRecord == null)primaryRecord = record;

						if (record.article.publicationDate == null)
							record.article.publicationDate = new Date();

						if (primaryRecord != record && primaryRecord.newIsBetter(record))
							primaryRecord = record;
					}

					for (ExperimentalProperty record : similarRecords) // list should include also our record
					{
						record.firstEntry = primaryRecord.id; // if only hidden exist -- reference to itself
						Globals.session().saveOrUpdate(record);
						if (recordIdSet.contains(record.id))
							ids.remove(record.id);
					}

				}

			}

			if (lastRestartSize - ids.size() >= 100)
			{
				lastRestartSize = ids.size();
				log("Primary record links processed for " + (totalSize - ids.size()) + " out of "+totalSize+" records in "+ ((int)((System.nanoTime() - timer)/1E7))/100. +" seconds");
				Globals.restartMainTransaction(true);
			}
		}

		Globals.commitMainTransaction();
	}

	@SuppressWarnings("unchecked")
	public List<ExperimentalProperty> getSimilarRecords(ExperimentalProperty ep)
	{

		Criteria duplicatesCriteria =  Globals.session().createCriteria(ExperimentalProperty.class);

		// not deleted
		duplicatesCriteria.add(Restrictions.isNull("deleted"));

		// Scope: include records, visible to the introducer of the current record
		//duplicatesCriteria.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE)); // only publicly available records

		// Property
		duplicatesCriteria.add(Restrictions.eq("property", ep.property));

		// Molecule
		duplicatesCriteria.createAlias("molecule", "mol");
		duplicatesCriteria.add(Restrictions.eq("mol.mapping1", ep.molecule.mapping1));

		// Value
		if (ep.property.isQualitative())
			duplicatesCriteria.add(Restrictions.eq("option", ep.option));
		else
		{
			duplicatesCriteria.add(Restrictions.or(Restrictions.eq("value", ep.value), Restrictions.eq("canonicalValue", ep.canonicalValue)));
		}

		duplicatesCriteria.addOrder(Order.asc("id"));

		return duplicatesCriteria.list();
	}
}


/*
@DatabaseMaintenanceJob
public class FirstEntryRecalculateTask extends OchemCronjobTask {

	public static void main(String[] args) 
	{
		new FirstEntryRecalculateTask().executeTask();
	}

	@SuppressWarnings("unchecked")
	private void resetReferencingRecords(ExperimentalProperty ep)
	{
		try 
		{
			List<ExperimentalProperty> referencingRecords = 
					Globals.session().createCriteria(ExperimentalProperty.class).add(Restrictions.or(
							Restrictions.eq("firstEntry", ep.id),
							Restrictions.eq("firstEntryIncludingHidden", ep.id)
							)).list();

			log(referencingRecords.size() + " referencing records found for record R" + ep.id);

			for (ExperimentalProperty referencingRecord : referencingRecords) 
			{
				referencingRecord.firstEntry = null;
				referencingRecord.firstEntryIncludingHidden = null;
				log("Invalidating a referenced record " + referencingRecord.id);
			}

		} catch (Exception e) {
			log("ERROR while fetching referencing records for ep id " + ep.id);
			e.printStackTrace();
		}
	}

	@SuppressWarnings("unchecked")
	public void executeTask()
	{	
		boolean FIX = true;

		if(FIX) return;

		ThreadScope.get().disableTrackChanges = true;
		long timer = System.nanoTime();
		Globals.startMainTransaction();

		List<Long> recordIds = 	Globals.session().createCriteria(ExperimentalProperty.class)
				.add(Restrictions.or(Restrictions.isNull("firstEntry"), Restrictions.isNull("firstEntryIncludingHidden")))
				.add(Restrictions.ne("molecule.id", 6933L)) // no empty molecules
				.createAlias("owner", "o")
				.createAlias("introducer", "i")
				.add(Restrictions.not(Restrictions.like("o.login", OCHEMTestHelpers.TEST_USER_PREFIX+"%")))
				.add(Restrictions.not(Restrictions.like("i.login", OCHEMTestHelpers.TEST_USER_PREFIX+"%")))
				.setProjection(Projections.id())
				.list();

		Set<Long> recordIdSet = new HashSet<Long>();
		recordIdSet.addAll(recordIds);

		int lastRestartSize, totalSize;
		lastRestartSize = totalSize = recordIds.size();

		log("Total number of records for primary record calculation: " + totalSize);

		while (recordIds.size() > 0)
		{
			Long recordId = recordIds.remove(0);
			ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, recordId);

			if (ep != null && ep.firstEntry != null && ep.firstEntryIncludingHidden != null) // ???
				continue;

			resetReferencingRecords(ep);

			//List<ExperimentalProperty> similarRecords = getSimilarRecords(ep); // IVT TODO -- disabled to stop infinite calculations

			List<ExperimentalProperty> similarRecords = new ArrayList<ExperimentalProperty>(); similarRecords.add(ep);

			ExperimentalProperty primaryRecord = ep;
			ExperimentalProperty primaryRecordIncludingHidden = ep;

			for (ExperimentalProperty record : similarRecords)
			{
				if (record.article.publicationDate == null)
					continue;

				if (record.rights == Globals.RIGHTS_FREELY_AVAILABLE)
					if (primaryRecord.article.publicationDate != null && !primaryRecord.isDeleted())
					{
						if (record.article.publicationDate.before(primaryRecord.article.publicationDate))
							primaryRecord = record;
						else if (record.article.publicationDate.equals(primaryRecord.article.publicationDate) && (record.id < primaryRecord.id) )
							primaryRecord = record;
					}
					else
						primaryRecord = record;

				if (primaryRecordIncludingHidden.article.publicationDate != null && !primaryRecordIncludingHidden.isDeleted()) // deleted records are always "less primary"
				{
					if (record.article.publicationDate.before(primaryRecordIncludingHidden.article.publicationDate))
						primaryRecordIncludingHidden = record;
					else if (record.article.publicationDate.equals(primaryRecordIncludingHidden.article.publicationDate) && (record.id < primaryRecordIncludingHidden.id) )
						primaryRecordIncludingHidden = record;
				}
				else
					primaryRecordIncludingHidden = record;
			}

			//Similar records include the record itself as well
			//Inserted piece to not count on it
			ep.firstEntry =  primaryRecord.id;
			ep.firstEntryIncludingHidden = primaryRecordIncludingHidden.id;
			Globals.session().saveOrUpdate(ep);
			log("Primary record for " + ep.id + " is set to " + primaryRecord.id + "/" + primaryRecordIncludingHidden.id);
			///
			for (ExperimentalProperty record : similarRecords)
			{
				record.firstEntry = primaryRecord.id;
				record.firstEntryIncludingHidden = primaryRecordIncludingHidden.id;
				Globals.session().saveOrUpdate(record);
				log("Primary record for " + record.id + " is set to " + primaryRecord.id + "/" + primaryRecordIncludingHidden.id);

				if (recordIdSet.contains(record.id))
					recordIds.remove(record.id);
			}

			if ( lastRestartSize - recordIds.size() >= 10 )
			{
				lastRestartSize = recordIds.size();
				log("Primary record links processed for " + (totalSize - recordIds.size()) + " out of "+totalSize+" records in "+ (System.nanoTime() - timer)/1E9 +" seconds");
				Globals.restartMainTransaction(true);
			}
		}

		Globals.commitMainTransaction();
	}

	@SuppressWarnings("unchecked")
	public List<ExperimentalProperty> getSimilarRecords(ExperimentalProperty ep)
	{
		//Criteria relatedRecordsCriteria = Globals.session().createCriteria(ExperimentalProperty.class);

		Criteria duplicatesCriteria =  Globals.session().createCriteria(ExperimentalProperty.class);

		// Scope: include records, visible to the introducer of the current record
		if (ep.introducer == null)
			duplicatesCriteria.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE));
		else
			ExperimentalProperty.addAccessRestrictions(duplicatesCriteria, Globals.RIGHTS_READ, ep.introducer, false, true);
		duplicatesCriteria.add(Restrictions.isNull("deleted"));

		// Property
		duplicatesCriteria.add(Restrictions.eq("property", ep.property));

		// Molecule
		if(ep.molecule != null)
		{
			duplicatesCriteria.createAlias("molecule", "mol");
			duplicatesCriteria.add(Restrictions.eq("mol.mapping1", ep.molecule.mapping1));
		}
		else
			duplicatesCriteria.add(Restrictions.isNull("molecule"));

		// Value
		if (ep.property.isQualitative())
			duplicatesCriteria.add(Restrictions.eq("option", ep.option));
		else
		{
			duplicatesCriteria.add(Restrictions.or(Restrictions.eq("value", ep.value), Restrictions.eq("canonicalValue", ep.canonicalValue)));
			//if(ep.unit != null)
			//	duplicatesCriteria.add(Restrictions.eq("unit", ep.unit));
		}

		// Find either dublicates, or referencing records
		//duplicatesCriteria.add(Restrictions.or(duplicatesCriteria, Restrictions.eq("firstEntry", ep.id)));
		duplicatesCriteria.addOrder(Order.asc("id"));

		return duplicatesCriteria.list();
	}
}*/
