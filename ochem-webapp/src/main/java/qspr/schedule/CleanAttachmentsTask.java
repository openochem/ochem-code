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

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.type.IntegerType;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.DataReferenceCleaner;
import qspr.metaserver.transport.NoSQLCleanerThread;
import qspr.metaserver.transport.NoSqlTransport;

@DatabaseMaintenanceJob
public class CleanAttachmentsTask extends OchemCronjobTask implements DataReferenceCleaner {

	private static final Logger logger = LogManager.getLogger(CleanArticleTask.class);

	@SuppressWarnings("unchecked")
	public void executeTask() throws Exception
	{	
		ThreadScope.get().disableTrackChanges = true;
		Globals.startMainTransaction();
		logger.info("Selecting valid attachment IDs");
		List<Integer> real_ids = Globals.session().createSQLQuery(
				"select teacher_task_template from ModelTemplate union distinct " +
						"select applier_task_template from ModelTemplate union distinct " +
						"select attachment_id from ArticleUserPdf union distinct " +
						"select attachment from Model union distinct " +
						"select microattachment from Model union distinct " +
						"select ready_model_attachment from Model union distinct " +
						"select calculated_descriptors from Model union distinct " +
						"select attachment_id from PendingTask union distinct " +
						"select stat_original from ModelMapping union distinct " +
						"select stat_recalculated from ModelMapping union distinct "+
						"select attachment_id from ShuffleKey union distinct "+
						"select file_attachment_id from BatchUpload union distinct " +
						"select ready_task from PendingTask union distinct " +
						"select attachment_id from UserAttachment union distinct " +
						"select configuration from ModelConfigurationTemplate" 
				).list();
		logger.info("Selecting all attachment IDs");
		List<Integer> all_ids = Globals.session().createSQLQuery("select attachment_id from Attachment").list();
		logger.info(String.format("Substracting two lists of with IDs %d and %d elements", all_ids.size(), real_ids.size()));

		Globals.commitMainTransaction();
		all_ids.removeAll(real_ids); //We now have only orphan attachment ids here
		Globals.startMainTransaction();

		logger.info("Found " + all_ids.size() + " orphan attachments");
		while (all_ids.size() > 0)
		{
			List<Integer> del_ids = all_ids.subList(0, Math.min(50, all_ids.size()));
			// Caution - "in" queries do not work well if there are a lot of IDs, that is a lot of orphan attachments (lets say thousands) // Midnighter on Oct 19, 2011
			logger.info("Deleting " + del_ids.size() + " attachments");
			Globals.session().createSQLQuery("delete from Attachment where attachment_id in (:list)").setParameterList("list", del_ids, IntegerType.INSTANCE).executeUpdate(); 
			logger.info("Deleted "+del_ids.size()+" orphan attachments: "+del_ids);
			del_ids.clear();
			Globals.restartMainTransaction(true);
		}
		log("No more orphan attachments found");
		Globals.commitMainTransaction();

		NoSQLCleanerThread clean = new NoSQLCleanerThread("Attachment", NoSqlTransport.DEFAULT_COLLECTION, this);
		clean.cleanCollection();
	}

	public static void main(String[] args)
	{
		OCHEMConfiguration.mirror = false;
		new CleanAttachmentsTask().executeInternal();
	}

	@Override
	public Set<String> getValidTaskReferenceIds() {
		Set<String> keepIds = new HashSet<String>();
		synchronized (this)
		{
			try
			{
				Globals.startMainTransaction();
				@SuppressWarnings("unchecked")
				List<String> attachmentIdsList = Globals.session().createSQLQuery("select attachment_md5 from Attachment").list();
				Globals.commitMainTransaction();

				for(String ref:attachmentIdsList) // we will not delete old reference style data
					if(ref != null && ref.length() == Task.REFERENCE_LENGTH) // otherwise it is an old style reference: will remain intact
						keepIds.add(ref);
				log("keeping "+keepIds.size());
			} catch (Exception e)
			{
				throw e;
			}
		}
		return keepIds;
	}

	@Override
	public void taskToBeDeleted(String collection, int size) {
		// TODO Auto-generated method stub
	}

	@Override
	public Logger getLogger() {
		return logger;
	}
}
