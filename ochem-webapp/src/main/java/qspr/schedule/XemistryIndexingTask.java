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
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.exception.ConstraintViolationException;
import org.hibernate.exception.DataException;
import org.hibernate.exception.GenericJDBCException;
import org.hibernate.exception.JDBCConnectionException;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.workflow.utils.QSPRConstants;
import qspr.util.WrapperThread;

/**
 * A cron job that adds the OCHEM molecules to the Xemistry (Cactvs) index
 * 
 * @author midnighter
 *
 */
@DatabaseMaintenanceJob
public class XemistryIndexingTask extends OchemCronjobTask
{
	private static final Logger logger = LogManager.getLogger(XemistryIndexingTask.class);

	static Integer start = 0 , maxsize = 0;

	static boolean debug = false;

	/**
	 * The list of "bad" molecule IDs that could not be indexed and will be skipped on the next runs
	 */
	private Set<Integer> badIDs = new TreeSet<Integer>();

	private String query ="insert ignore into StructureQuery select\n " + 
			"	mapping2_id,\n" + 
			"	ens_string_property(uncompress(molecule_data),'E_QUERY_SCREEN',null,'rawbinary'),\n" +
			"	ens_blob_property(uncompress(molecule_data),'E_SCREEN'),\n" +
			"	ens_string_property(uncompress(molecule_data),'E_MINIMOL',null,'rawbinary')\n" + 
			"	from Molecule left join StructureQuery using (mapping2_id) where StructureQuery.mapping2_id is null and molecule_id in (:ids) ";

	public void executeTask() throws Exception 
	{

		new WrapperThread()
		{
			@Override
			@SuppressWarnings("unchecked")
			public void wrapped() throws Exception
			{
				List<Object[]> rows = Globals.session().createSQLQuery("select max(molecule_id) from Molecule")
						.list();

				maxsize = Integer.parseInt(""+rows.get(0));

				boolean first  = start == 0;

				if(maxsize > QSPRConstants.INHOUSE_XEMISTRY){
					if(first){
						rows = Globals.session().createSQLQuery("select max(molecule_id) from Molecule where molecule_id < " + QSPRConstants.INHOUSE_XEMISTRY)
								.list();
						maxsize = Integer.parseInt(""+rows.get(0)); // up to maximum number of not inhouse values
					}else
						start = start > QSPRConstants.INHOUSE_XEMISTRY ? start : QSPRConstants.INHOUSE_XEMISTRY ; // Reset start to work with inhouse values
				}

				List<Integer> unindexedIds = getUnindexedIDs();
				while (unindexedIds.size() > 0)
				{
					try
					{
						indexMolsTolerantly(unindexedIds);
						Globals.restartAllTransactions(true);
						unindexedIds = getUnindexedIDs();
					}
					catch (JDBCConnectionException e)
					{
						logger.warn("JDBC Exception occured. Retrying.");
						throw e;
					}
				}

				if(first && badIDs.size() >0)
					logger.error("Number of failed molecules " + badIDs.size() + " >= XEMISTRY_INDEXING_TASK = " + 
							QSPRConstants.XEMISTRY_INDEXING_TASK + ". No idexing will be performed. Failed entries are:\n\n" + badIDs);

			}
		}.run();
	}

	@Override
	protected boolean shouldRun()
	{
		if (!OCHEMConfiguration.xemistryIndexing) return false;

		return super.shouldRun();
	}

	@SuppressWarnings("unchecked")
	protected List<Integer> getUnindexedIDs()
	{
		logger.info("Quering unindexed molecules...");

		List<Integer> list = new ArrayList<Integer>();
		Set<Integer> unindexedIds = new HashSet<Integer>();

		for(int i = 0; start < maxsize ; start += QSPRConstants.XEMISTRY_INDEXING_TASK, i++){

			if(i % 10 == 0){
				Globals.restartAllTransactions(false);  // to avoid timeout
				logger.info(" Analysing\t" + start + " ouf of " + maxsize + " processed " + (int)(10000.*start/maxsize)/100.+"% size:" + unindexedIds.size());
			}

			if(debug)logger.info("Analysing\t" + start + " ouf of " + maxsize);

			//			String query = "select mapping2_id from ExperimentalProperty natural join Molecule left join StructureQuery " + 
			//					"using (mapping2_id) where StructureQuery.mapping2_id is null limit " + start + ", " + QSPRConstants.XEMISTRY_INDEXING_TASK;

			String query = "select distinct molecule_id from Molecule natural join ExperimentalProperty left join StructureQuery using(mapping2_id) where molecule_id >= " + 
					start + "\n and molecule_id < "+  (start + QSPRConstants.XEMISTRY_INDEXING_TASK) + " and StructureQuery.mapping2_id is NULL";

			if(debug)logger.info(query);

			List<Integer> ids = Globals.session().createSQLQuery(query).list();

			if(ids.size() == 0)continue;

			unindexedIds.addAll(ids);
			if(unindexedIds.size() > QSPRConstants.XEMISTRY_INDEXING_TASK)
				unindexedIds.removeAll(badIDs);

			if(unindexedIds.size() > QSPRConstants.XEMISTRY_INDEXING_TASK){
				list.addAll(unindexedIds);
				logger.info("Will index " + unindexedIds.size() + " start=" + start);
				return list;
			}
		}

		start = maxsize;
		list.addAll(unindexedIds);
		logger.info("Will finally index " + unindexedIds.size());
		return list;

	}

	/**
	 * Smart logic of indexing: calculate in batches, switch to one-by-one processing in case of failures
	 * @param ids
	 */
	private void indexMolsTolerantly(List<Integer> ids)
	{
		logger.info(String.format("Indexing %d molecules (%d skipped)", ids.size(), badIDs.size()));
		if (!indexMolsSafely(ids))
		{
			Globals.restartAllTransactions(false);
			for (Integer id : ids)
			{
				if (!indexMolsSafely(Arrays.asList(new Integer[]{id})) && badIDs.isEmpty())
					logger.warn(" Molecule M" + id + " skipped. Other errors will not be reported. Query was:\n" + query + "\n");
				badIDs.add(id);
			}
		}
	}

	/**
	 * Index a list of molecules
	 * Returns true on success and false otherwise 
	 * @param ids
	 * @return
	 */
	private boolean indexMolsSafely(List<Integer> ids)
	{
		try
		{
			indexMols(ids);
			return true;
		}
		catch (GenericJDBCException e)
		{
			return false;
		}
		catch (DataException e)
		{
			return false;
		}
		catch (ConstraintViolationException e)
		{
			return false;
		}
	}

	/**
	 * Add the specified molecules to the cactvs index
	 * @param ids - identifiers of the molecules to add
	 */
	protected void indexMols(List<Integer> ids)
	{
		// As provided by the Wolf Ihlenfeldt
		Globals.session().createSQLQuery(query)
		.setParameterList("ids", ids)
		.executeUpdate();
	}


	public static void main(String[] args) throws Exception
	{
		XemistryIndexingTask t = new XemistryIndexingTask();
		//debug = true;
		t.executeTask();
	}


}
