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

package qspr.metaserver.transport;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import qspr.util.DaemonThread;

import com.eadmet.utils.mailer.Mailer;

public class NoSQLCleanerThread extends DaemonThread
{
	private static final long WAIT_CLEANING_MONGODB = 15 * 60 * 1000;
	private String database;
	private String collection;

	DataReferenceCleaner cleaner;

	public NoSQLCleanerThread(String database, String collection, DataReferenceCleaner cleaner)
	{
		super(cleaner.getLogger());
		this.database = database;
		this.collection = collection;
		this.cleaner = cleaner;
	}

	@Override
	public void wrapped() throws InterruptedException
	{
		try
		{
			Thread.sleep(WAIT_CLEANING_MONGODB);
			cleanCollection();
		} catch (Exception e)
		{
			log.info(e);
			Mailer.notifyDevelopers(e, "Metaserver$MaintenanceThread/NoSQLclean");
		}
	}

	public void cleanCollection() throws Exception
	{
		log.info("[Clean MongoDB] Quering MongoDB IDs for " + database+"/"+collection);
		Set<String> nosqlIds = NoSqlTransport.getReferenceIds(database, collection);
		log.info("[Clean MongoDB] Query for " + database+"/"+collection + " completed. Now waiting for " + 
				WAIT_CLEANING_MONGODB / (60* 1000) + "  minutes before performing clean-up...");

		Thread.sleep(WAIT_CLEANING_MONGODB);

		log.info("[Clean MongoDB] Quering valid IDs for " + collection);
		Set<String> databaseIds = cleaner.getValidTaskReferenceIds();

		log.info("[Clean MongoDB] Substracting two lists of identifiers of sizes " + nosqlIds.size() + " and " + databaseIds.size() + " for " + database+"/"+collection + "...");
		nosqlIds.removeAll(databaseIds);
		int size = nosqlIds.size();
		cleaner.taskToBeDeleted(collection, size);
		log.info("[Clean MongoDB] To delete in collection " + database+"/"+collection + ": " + size);
		long timer = System.nanoTime();

		List<String> setList=new ArrayList<String>(nosqlIds); 

		while (setList.size() > 0)
		{
			List<String> oneBatchIDs = setList.subList(0, Math.min(100, setList.size()));
			log.info("[Clean MongoDB] Deleting " + oneBatchIDs.size() + " records out of " + setList.size() + " in collection " + database+"/"+collection);
			NoSqlTransport.deleteReferenceByIds(database, collection, oneBatchIDs);
			oneBatchIDs.clear();
			cleaner.taskToBeDeleted(collection, setList.size());
		}
		timer = (System.nanoTime() - timer) / (1000 * 1000 * 1000);
		log.info("[Clean MongoDB] Deletion of outdated MongoDB records for " + database+"/"+collection + " completed, took " + timer + " seconds for " + size + " records.");
	}
}