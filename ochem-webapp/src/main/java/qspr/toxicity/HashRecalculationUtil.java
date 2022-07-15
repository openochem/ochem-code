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

package qspr.toxicity;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.entities.ExperimentalProperty;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.mailer.Mailer;


public class HashRecalculationUtil
{
	public static void emptyHashRecalculationFull()
	{
		new EmptyHashRecalculatorFull().doIterate();
	}

	public static void allHashRecalculationLightweight()
	{
		new AllHashRecalculatorLightweight().doIterate();
	}

	public static void allHashRecalculation()
	{
		ThreadScope.get().disableTrackChanges = true;
		Globals.startMainTransaction();
		Globals.session().createQuery("update ExperimentalProperty set md5 = null").executeUpdate();
		Globals.session().createQuery("update ExperimentalProperty set ep_status = null where ep_status = 0").executeUpdate();
		Globals.session().flush();
		Globals.commitMainTransaction();
	}


	public static void main(String... args)
	{
		//new EmptyCanonicalValueCalculator().doIterate();
		new AllHashRecalculatorLightweight().doIterate();
	}
}

abstract class EPHashCalculator
{
	private Logger logger = LogManager.getLogger(getClass());

	private static final int BATCHSIZE = 100;

	EpHashCalculationTemplates hct = new EpHashCalculationTemplates();
	EpCriteriaTemplates ect = new EpCriteriaTemplates();


	public boolean safeRun = false;

	@SuppressWarnings({ "deprecation", "unchecked" })
	List<ExperimentalProperty> getAllLazyFromSublit(List<Long> sublist){
		// lazy loading everything required for the analyzed records
		Criteria epc = Globals.session().createCriteria(ExperimentalProperty.class)
				.createAlias("molecule", "mol", Criteria.LEFT_JOIN)
				.createAlias("mol.mapping1", "mp1", Criteria.LEFT_JOIN)
				.createAlias("mol.mapping2", "mp2", Criteria.LEFT_JOIN)
				.createAlias("moleculenames", "mn", Criteria.LEFT_JOIN)
				.createAlias("property", "prop", Criteria.LEFT_JOIN)
				.createAlias("conditions", "cond", Criteria.LEFT_JOIN)
				.add(Restrictions.in("id", sublist))
				.setResultTransformer(Criteria.DISTINCT_ROOT_ENTITY);

		return epc.list();
	}

	@SuppressWarnings("unchecked")
	public void doIterate()
	{
		ThreadScope.get().disableTrackChanges = true;
		Globals.startMainTransaction();
		List<Long> ids = new ArrayList<Long>();
		String name = "" + this; name = name.substring(0,name.indexOf('@')); name = name.substring(name.indexOf('.')+1);name = name.substring(name.indexOf('.')+1);

		//ids.add(27982761l);

		Criteria c;
		c = getRecordsCriteria()
				//.add(Restrictions.ge("id",27982761l))
				.add(Restrictions.isNull("deleted"));

		if(ids.size() == 0)
			ids = c.list(); // if no pre-set list

		Collections.sort(ids);

		int size = ids.size();

		if(size>0)
			logger.info("Got Total " + ids.size() + " initial records to update in range [" + ids.get(0) + " ; " + ids.get(size-1) + "] "  + name);

		long global = Calendar.getInstance().getTimeInMillis(), timeBatch = 0;
		safeRun = false; //With manual checks and so on
		int  failedEntriesToProcess = 0, attempt = 0;

		while (ids.size() > 0)
		{
			long time = Calendar.getInstance().getTimeInMillis();
			int max = safeRun ? 1 : BATCHSIZE; 
			List<Long> sublist = ids.subList(0, Math.min(max,ids.size()));

			try
			{
				List<ExperimentalProperty> eps = getAllLazyFromSublit(sublist);

				if(attempt >= 1) { // attempt to manually setting md5 to the second record to make possible indexing and add second record
					eps = interchangeTwoFailing(eps.get(0));
					logger.info("manually changing orders for records: " + eps.get(0) + " and "  + eps.get(1) + " " + name);

					if(attempt == 2) { // does not work; the order should be restored and records re-fetched
						ExperimentalProperty ep = eps.get(0);
						ep.md5 = null;
						Globals.session().saveOrUpdate(ep);							
						Globals.restartMainTransaction(true);
						eps = getAllLazyFromSublit(Arrays.asList(ep.id,eps.get(1).id));
						eps.add(0, eps.get(1));
						eps.remove(2);
					}
				}


				for (ExperimentalProperty ep : eps)
					doAction(ep);

				Globals.restartMainTransaction(true);
				failedEntriesToProcess -= sublist.size();
				sublist.clear();

			} catch (Exception e)
			{
				Globals.restartMainTransaction(false);
				logger.info("Restarting batch in safe mode due to exception: " + e.getMessage() + " " + name);
				if (safeRun) {
					logger.info(e);					
					if(attempt++ <= 2) continue;
					System.out.println("Still does not work");
					Mailer.notifyDevelopers("Hash calculation failed in safeRun", "Hash recalculation failed for R" +sublist.get(0) + " with exception\n"+OCHEMUtils.exceptionToString(e));
					return;
				}

				attempt = 0;
				safeRun = true;
				timeBatch = Calendar.getInstance().getTimeInMillis();
				failedEntriesToProcess = sublist.size();
				continue;
			}

			attempt = 0;
			sublist.clear();

			long delta = (Calendar.getInstance().getTimeInMillis() - time)/1000; // in sec
			long globaldelta = (Calendar.getInstance().getTimeInMillis() - global)/1000; // in sec

			String message = String.format("Processed %d out of %d, total time: %ds, last batch: %ds, remaining: %ds", 
					size - ids.size(), size, globaldelta, delta, globaldelta * ids.size()/(size - ids.size()));

			logger.info(message + " " + name);

			if(safeRun && failedEntriesToProcess <= 0 && ids.size() > 0) {
				delta =  (Calendar.getInstance().getTimeInMillis() - timeBatch)/1000;
				logger.info("Time to execute batch was " +  delta + "s. Restoring batch size to the default value: " + BATCHSIZE + " " + name);
				logger.info("Still have " + ids.size() + " records to update in range [" + ids.get(0) + " ; " + ids.get(ids.size()-1) + "] "  + name);
			}

			safeRun = failedEntriesToProcess > 0;
		}

		Globals.commitMainTransaction();

	}

	protected List<ExperimentalProperty> interchangeTwoFailing(ExperimentalProperty ep){
		ep.ep_status = null; // cleaning the status, otherwise duplicates will not be found
		ep.updateHash();
		Globals.restartMainTransaction(false);
		List<ExperimentalProperty> duplicates = ep.getAllDublicates();
		ExperimentalProperty duplicate = duplicates.get(0);
		duplicate.md5 = null;
		Globals.session().saveOrUpdate(duplicate);
		Globals.restartMainTransaction(true);
		List<Long> ids2 = new ArrayList<Long>();
		ids2.add(ep.id);
		ids2.add(duplicate.id);
		return getAllLazyFromSublit(ids2);
	}

	public abstract Criteria getRecordsCriteria();
	public abstract void doAction(ExperimentalProperty ep);

}

class MixtureHashRecalculatorFull extends EPHashCalculator 
{

	@Override
	public Criteria getRecordsCriteria() 
	{
		return ect.getMixtureRecords();
	}

	@Override
	public void doAction(ExperimentalProperty ep) 
	{
		if (safeRun)
			hct.calculateHashFullSlow(ep);
		else
			hct.calculateHashFullQuick(ep);

	}

}


class AllHashRecalculatorFull extends EPHashCalculator 
{

	@Override
	public Criteria getRecordsCriteria() 
	{
		return ect.getAllValidRecords();
	}

	@Override
	public void doAction(ExperimentalProperty ep) 
	{
		if (safeRun)
			hct.calculateHashFullSlow(ep);
		else
			hct.calculateHashFullQuick(ep);

	}

}

class EmptyHashRecalculatorFull extends EPHashCalculator 
{

	@Override
	public Criteria getRecordsCriteria() 
	{
		return ect.getEmptyHashValidRecords();
	}

	@Override
	public void doAction(ExperimentalProperty ep) 
	{
		if (safeRun)
			hct.calculateHashFullSlow(ep);
		else
			hct.calculateHashFullQuick(ep);

	}

}

class EmptyCanonicalValueCalculator extends EPHashCalculator 
{
	@Override
	public Criteria getRecordsCriteria() 
	{
		return ect.getSomeEmptyCanonicalValueRecords();
	}

	@Override
	public void doAction(ExperimentalProperty ep) 
	{
		hct.calculateHashLightweight(ep);		
	}

}

class AllHashRecalculatorLightweight extends EPHashCalculator 
{
	@Override
	public Criteria getRecordsCriteria() 
	{
		return ect.getAllValidRecords();
	}

	@Override
	public void doAction(ExperimentalProperty ep) 
	{
		if (!safeRun)
			hct.calculateHashLightweight(ep);
		else
			hct.calculateHashFullSlow(ep);
	}

}


class EpCriteriaTemplates
{
	public Criteria getEmptyHashValidRecords() // Records that should have a hash (valid and not deleted) but don't
	{
		return 	Globals.session().createCriteria(ExperimentalProperty.class)
				.add(Restrictions.isNull("md5"))
				.add(Restrictions.isNull("deleted"))
				.add(Restrictions.or(
						Restrictions.and(
								Restrictions.ne("ep_status", Integer.valueOf(ExperimentalProperty.STATUS_ERROR)), 
								Restrictions.ne("ep_status", Integer.valueOf(ExperimentalProperty.STATUS_INVALID))), 
						Restrictions.isNull("ep_status")))
				.setProjection(Projections.id());
	}

	public Criteria getMixtureRecords() // one time mixture runthrough... can be changed to suit the other needs later
	{
		return 	Globals.session().createCriteria(ExperimentalProperty.class)
				.createAlias("conditions", "c")
				.createAlias("c.values", "v")
				.createAlias("v.property", "p")
				.add(Restrictions.eq("p.name", QSPRConstants.MIXTURE_CONDITION))
				.setProjection(Projections.id());
	}

	public Criteria getSomeEmptyCanonicalValueRecords() 
	{
		return 	Globals.session().createCriteria(ExperimentalProperty.class)
				.add(Restrictions.isNull("deleted"))
				.add(Restrictions.gt("value", 0D))
				.add(Restrictions.eq("canonicalValue", 0D))
				.createAlias("property", "p")
				//.add(Restrictions.like("p.name","Retention Factor"))
				.setMaxResults(10)
				.setProjection(Projections.id());
	}

	public Criteria getAllValidRecords() // Records that should have a hash (valid and not deleted) but don't
	{
		return 	Globals.session().createCriteria(ExperimentalProperty.class)
				.add(Restrictions.isNull("deleted"))
				.setProjection(Projections.id());
	}
}

class EpHashCalculationTemplates
{
	private static transient final Logger logger = LogManager.getLogger();

	public void calculateHashFullQuick(ExperimentalProperty ep)
	{
		try
		{
			ep.updateHash();
		}
		catch (Exception e)
		{
			ep.ep_status = ExperimentalProperty.STATUS_INVALID;
			ep.errorComment = e.getMessage();
			logger.info("(ERR) R"+ep.id+": "+ep.errorComment);
			ep.updateHash();
		}
		Globals.session().saveOrUpdate(ep);
	}

	public void calculateHashFullSlow(ExperimentalProperty ep)
	{
		ep.ep_status = null; // cleaning the status, otherwise duplicates will not be found
		calculateHashFullQuick(ep);
		ep.resolveConflictsAndSave();
		Globals.restartMainTransaction(true);
		if (ep.errorComment != null)
			logger.info("(DBL) R"+ep.id+": "+ep.errorComment);
		Globals.session().flush();
	}

	private static boolean equal(String s1, String s2)
	{
		if (s1 == null)
			return (s2 == null);
		else
			return s1.equals(s2);
	}

	public void calculateHashLightweight(ExperimentalProperty ep)
	{
		String oldMD5 = ep.hash;
		try
		{
			ep.ep_status = null;
			ep.updateHash();
		}
		catch (Exception e)
		{
			ep.ep_status = ExperimentalProperty.STATUS_INVALID;
			ep.errorComment = e.getMessage();
			ep.updateHash();
		}

		if (equal(oldMD5, ep.hash)) // Hash recalculation brought no change whatsoever... could be status or error message change, but it's unimportant
		{
			Globals.session().saveOrUpdate(ep);
			return;
		}

		ep.resolveConflictsAndSave();
		Globals.session().flush();

		if (equal(oldMD5, ep.md5)) // Hash recalculation brought some change, but it was negated by duplicate resolve check
		{
			Globals.session().saveOrUpdate(ep);
			Globals.session().flush();

			return;
		}

		logger.info("(CHG) "+oldMD5+"->"+ep.hash+ " for: " + ep);
		if (ep.errorComment != null)
			logger.info("(ERR) "+ep.errorComment);

		Globals.session().saveOrUpdate(ep);
		Globals.session().flush();
	}
}
