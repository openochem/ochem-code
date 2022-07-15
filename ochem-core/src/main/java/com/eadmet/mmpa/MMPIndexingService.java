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

package com.eadmet.mmpa;

import java.io.IOException;
import java.net.MalformedURLException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.SQLQuery;
import org.hibernate.criterion.Restrictions;
import org.hibernate.type.IntegerType;
import org.hibernate.type.LongType;
import org.hibernate.type.StringType;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.Mapping2;
import qspr.entities.Mapping2.MMPAIndex;
import qspr.entities.Molecule;
import qspr.metaserver.CalculationClient;
import qspr.metaserver.configurations.MMPFragConfiguration;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.workflow.utils.QSPRConstants
;
import qspr.util.WrapperThread;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.mmpa.domain.MMPFragment;
import com.eadmet.mmpa.domain.MMPIndex;
import com.eadmet.mmpa.domain.MMPTransformation;
import com.eadmet.mmpa.domain.MMPair;
import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.config.ConfigurableClass;
import com.eadmet.utils.config.ConfigurableProperty;

/**
 * MMPs require a computationally- and memory-intensive indexing of molecules.
 * 
 * This service implements indexing in two steps: 
 *  - creation of a fragment index and 
 *  - identification of matched pairs
 *  
 * @author midnighter
 *
 */
@ConfigurableClass(name = "mmp", comment = "Molecular matched pairs module")
@SuppressWarnings("unchecked")
public class MMPIndexingService
{

	@ConfigurableProperty(name = "enabled", comment = "If enabled, the system will run the molecular matched pairs indexing tasks")
	public static boolean enabled = false;

	private static final Logger logger = LogManager.getLogger(MMPIndexingService.class);

	private final static int PENDING_TASK_LIMIT = 200;
	private final static int BATCH_SIZE = 500;
	private static final long INFREQUENT_TRANSFORMATIONS = 0; // all transformations will be kept

	private static final String MMPAService = "MMPAService";

	public int reSubmitIndexingTasks(Basket b) throws IOException, ClassNotFoundException, InterruptedException {

		int i= 0;
		List<Integer> mmp2 = new ArrayList<Integer>();
		do{			
			List<Molecule> set = new ArrayList<Molecule>();
			for(; set.size() < BATCH_SIZE && i < b.entries.size(); i++){
				Molecule mol = b.entries.get(i).ep.molecule;
				if(mmp2.contains(mol.mapping2.id))continue;
				if(mol.mapping2.mmpaIndexStatus != MMPAIndex.MMP_INDEXED){ // unless already indexed
					set.add(mol);
					mmp2.add(mol.mapping2.id);
				}
			}
			if(set.size() != 0) submitIndexingTaskMolecules(set, TaskPriority.LOW);
			Globals.restartAllTransactions(true);
		}while(i < b.entries.size());

		return mmp2.size();
	}

	public void submitIndexingTasks() throws IOException, ClassNotFoundException, InterruptedException {

		int pendingFragTasks = getPendingTasksCount();
		logger.info("" + pendingFragTasks + " pending fragmentation tasks. The limit is 200.");

		List<Mapping2> mols;

		do{
			Globals.restartAllTransactions(true);
			mols = Globals.session().createCriteria(Mapping2.class)
					.add(Restrictions.eq("mmpaIndexStatus", MMPAIndex.MMP_SCHEDULED))
					.setMaxResults(BATCH_SIZE)
					.list();
			submitIndexingTask(mols, TaskPriority.LOW);
		}while(mols.size() > 0 && getPendingTasksCount() < PENDING_TASK_LIMIT);

	}

	private int submitIndexingTaskMolecules(List<Molecule> mols, int priority) throws IOException, ClassNotFoundException, InterruptedException {
		if(mols.size() == 0) return -1;

		DataTable dtMols = new DataTable();
		dtMols.addColumn(QSPRConstants.SDF_COLUMN);
		logger.info("Found " + mols.size() + "+ molecules for MMP fragmentation");
		for (Molecule mol : mols)
		{	
			AbstractDataRow r = dtMols.addRow();
			dtMols.setValue(0, mol.getData());
			r.addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, mol.mapping2.id);
		}

		for (Molecule mol : mols)
		{
			mol.mapping2.mmpaIndexStatus = MMPAIndex.MMP_SUBMITTED;
			Globals.session().saveOrUpdate(mol.mapping2);
		}

		logger.info("Posting " + dtMols.getRowsSize() + " molecules for MMP fragmentation");

		return new CalculationClient(MMPAService).postTask(new Task(QSPRConstants.MMPFrag, new MMPFragConfiguration(), new WorkflowNodeData(dtMols)).setPriority(priority).setMinRequiredMemory(4096));
	}

	private int submitIndexingTask(List<Mapping2> mp2s, int priority) throws IOException, ClassNotFoundException, InterruptedException {
		if(mp2s.size() == 0) return -1;

		DataTable dtMols = new DataTable();
		dtMols.addColumn(QSPRConstants.SDF_COLUMN);
		logger.info("Found " + mp2s.size() + "+ molecules for MMP fragmentation");
		for (Mapping2 mp2 : mp2s)
		{	
			AbstractDataRow r = dtMols.addRow();
			dtMols.setValue(0, mp2.getMolecule().getData());
			r.addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, mp2.id);
		}

		for (Mapping2 mp2 : mp2s)
		{
			mp2.mmpaIndexStatus = MMPAIndex.MMP_SUBMITTED;
			Globals.session().saveOrUpdate(mp2);
		}

		logger.info("Posting " + dtMols.getRowsSize() + " molecules for MMP fragmentation");

		return new CalculationClient(MMPAService).postTask(new Task(QSPRConstants.MMPFrag, new MMPFragConfiguration(), new WorkflowNodeData(dtMols)).setPriority(priority).setMinRequiredMemory(4096));
	}

	public void rebuildMMPs() {
		Globals.session().createQuery("delete from MMPair").executeUpdate();
		Globals.session().createQuery("delete from MMPTransformation").executeUpdate();

		List<MMPIndex> index = Globals.session().createCriteria(MMPIndex.class).list();
		logger.info("" + index.size() + " index elements found");
		int cnt = 0;
		for (MMPIndex mmpIndex : index)
		{
			processIndexElement(mmpIndex, true);
			if (cnt++ % 100 == 0)
			{
				logger.info("" + cnt + " elements processed");
				Globals.restartAllTransactions(true);
			}
		}
	}

	private int processIndexElement(MMPIndex index, boolean avoidDublicates) {

		int mmpCount = 0;
		String scaffold = index.scaffoldId;
		Long fragId = index.fragmentId;
		Integer mapping2 = index.mapping2Id;
		// Enumerate the matched pairs
		List<MMPIndex> indexElements = MMPIndex.getEntries(scaffold);

		for (MMPIndex entry : indexElements)
			if (entry.fragmentId != fragId)
			{
				if (avoidDublicates && entry.mapping2Id > mapping2)
					continue;
				if (mapping2.equals(entry.mapping2Id))
					continue;

				int[] molIDs = new int[]{mapping2, entry.mapping2Id};
				long[] fragIDs = new long[]{fragId, entry.fragmentId};
				Mapping2[] mp2s = new Mapping2[]{Repository.molecule.getMapping2(molIDs[0]), Repository.molecule.getMapping2(molIDs[1])};

				if (mp2s[0].mapping1.id.equals(mp2s[1].mapping1.id))
					continue;

				MMPair pair = MMPTransformation.getMMPair(fragIDs , mp2s);

				if (getForbiddenTransformationsIds().contains(pair.transformation.id))
					continue;

				pair.transformation.pairsCount++;
				Globals.session().save(pair);
				mmpCount++;
			}

		return mmpCount;
	}

	public int processTask(Task task) throws Exception {

		CalculationClient client = new CalculationClient(MMPAService);
		task.check();
		long time = Calendar.getInstance().getTimeInMillis();

		WorkflowNodeData wnd = WorkflowNodeData.fromTask(task);
		DataTable dtResult = wnd.ports.get(0);

		logger.info("Got " + dtResult.getRowsSize() + " index entries from MMP fragmentation. Updating the index..");

		long txTime = Calendar.getInstance().getTimeInMillis();
		int newEntries = 0, newMMPs = 0;

		for (dtResult.reset();dtResult.nextRow();)
		{
			Integer mapping2 = (Integer) dtResult.getCurrentRow().getAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM);
			Mapping2 mp2 = Repository.molecule.getMapping2(mapping2);

			if (dtResult.getCurrentRow().isError()) {
				mp2.mmpaIndexStatus = MMPAIndex.MMP_FAILED;
				continue;
			}

			for(int i=0;i*3 < dtResult.getColumnsSize();i++) {
				String scaffold = (String) dtResult.getValue(MMPFragConfiguration.SCAFFOLD_INCHI+i);
				if(scaffold == null)break; // no more data, stop
				String fragmentInchi = (String) dtResult.getValue(MMPFragConfiguration.FRAGMENT_INCHI+i);
				String fragmentSmiles = (String) dtResult.getValue(MMPFragConfiguration.FRAGMENT_SMILES+i);

				MMPFragment frag = MMPFragment.getFragment(fragmentInchi, fragmentSmiles);
				MMPIndex index = MMPIndex.saveIndex(scaffold, frag.id, mapping2);
				if (index != null) {
					newMMPs += processIndexElement(index, false);
					newEntries++;
				}
			}

			mp2.mmpaIndexStatus = MMPAIndex.MMP_INDEXED; // molecules was processed

			long interval = Calendar.getInstance().getTimeInMillis() - txTime;
			if (interval > 1000 * 5)
			{
				Globals.restartAllTransactions(true);
				if (newMMPs != 0)
					logger.info("Restarted transaction. " + "" + newMMPs + " new MMPs created in " + (Calendar.getInstance().getTimeInMillis() - time) + "ms. (" + ((Calendar.getInstance().getTimeInMillis() - time) / newMMPs) +" per MMP)");
				logger.info(MemoryUtils.memorySummary());
				txTime = Calendar.getInstance().getTimeInMillis();
			}
		}

		time = Calendar.getInstance().getTimeInMillis() - time;
		logger.info("" + newEntries + " new index entries created");
		if (newMMPs != 0){
			logger.info("" + newMMPs + " new MMPs created in " + time + "ms. (" + (time / newMMPs) +" per MMP)");
			MMPQueryService.startCacheUpdate(); /// is required to upload new indexes to cache
		}

		client.deleteTask(task.id);

		Globals.restartAllTransactions(true);

		return newEntries;
	}

	public void fetchIndexingTasks() throws Exception {
		CalculationClient client = new CalculationClient(MMPAService);
		Task task = client.getReadyTask(QSPRConstants.MMPFrag);

		int newEntriesTotal = 0;

		while (task != null) 
		{
			newEntriesTotal += processTask(task);
			task = client.getReadyTask(QSPRConstants.MMPFrag);
			Globals.restartAllTransactions(true);
		}

		if (newEntriesTotal > 0) {
			recalculatePairCounts();
			deleteInfrequentTransformations();
			deleteInvalidPairs();
		}
	}

	public void instantIndex(List<Mapping2> mp2s) throws Exception {

		ThreadScope.setStatus("Fragmentation and indexing of the molecule...", logger);
		Iterator<Mapping2> iter = mp2s.iterator();
		while (iter.hasNext())
			if (iter.next().mmpaIndexStatus == MMPAIndex.MMP_INDEXED)
				iter.remove();

		if (mp2s.isEmpty())
			return;

		CalculationClient client = new CalculationClient(MMPAService);
		int taskId = submitIndexingTask(mp2s, TaskPriority.EXTRA_HIGH);
		Task task = null;
		do {
			task = client.getTask(taskId);
			Thread.sleep(1000);
		} while (task == null);

		processTask(task);
	}

	public void instantIndex(Mapping2 molecule) throws Exception {
		List<Mapping2> list = new ArrayList<Mapping2>();
		list.add(molecule);
		instantIndex(list);
	}

	public MMPStatus getStatus(MMPQuery mmpQuery) {
		MMPStatus status = new MMPStatus();
		String[] idxStatusStr = new String[]{"Ignored", "Scheduled for indexation", "Submitted for fragmentation", "Indexed", "Failed"};
		SQLQuery query = null;
		if (mmpQuery == null)
			query = Globals.session().createSQLQuery("select mmpaIndexStatus, count(*) c from Mapping2 group by mmpaIndexStatus");
		else if (mmpQuery.basketId != null)
		{
			query = Globals.session().createSQLQuery("select mmpaIndexStatus, count(distinct mp2.mapping2_id) c from BasketEntry join ExperimentalProperty using (exp_property_id) join Molecule using (molecule_id) join Mapping2 mp2 using (mapping2_id) where (basket_id=:basketId) group by mmpaIndexStatus");
			query.setParameter("basketId", mmpQuery.basketId);
		}
		else if (mmpQuery.tagId != null)
		{
			query = Globals.session().createSQLQuery("select mmpaIndexStatus, count(distinct mp2.mapping2_id) c from MoleculeTag natural join Mapping1 natural join Mapping2 mp2 where (tag_id=:tagId) group by mmpaIndexStatus");
			query.setParameter("tagId", mmpQuery.tagId);
		}
		else
			throw new UserFriendlyException("Unsupported MMP query");

		List<Object[]> rows = query
				.addScalar("mmpaIndexStatus", IntegerType.INSTANCE)
				.addScalar("c", LongType.INSTANCE)
				.list();

		for (Object[] row : rows)
			status.molCount.put(idxStatusStr[((Integer)row[0])], (Long) row[1]);

		if (mmpQuery == null) {
			String[] tables = new String[]{"MMPTransformation", "MMPFragment", "MMPIndex", "MMPair"};
			for (String table : tables)
				status.getTableInfo(table).count = getRowsCount(table);

			rows = Globals.session().createSQLQuery("show table status")
					.addScalar("Name", StringType.INSTANCE)
					.addScalar("Data_length", LongType.INSTANCE)
					.addScalar("Index_length", LongType.INSTANCE)
					.list();
			for (Object[] row : rows)
				for (String table : tables)
					if (row[0].equals(table))
					{
						TableInfo info = status.getTableInfo(table);
						info.dataSize = (Long)row[1] / (1024 * 1024);
						info.indexSize = (Long)row[2] / (1024 * 1024);
						status.totalSize += info.size = info.dataSize + info.indexSize;
					}
		}

		return status;
	}

	public void scheduleForIndexation(MMPQuery query) {
		if (query.tagId != null)
			Globals.session().createSQLQuery("update Mapping2, MoleculeTag set mmpaIndexStatus=:mmpnew where mmpaIndexStatus=:mmpold and Mapping2.mapping1_id=MoleculeTag.mapping1_id and tag_id=" + query.tagId)
			.setParameter("mmpnew", MMPAIndex.MMP_SCHEDULED)
			.setParameter("mmpold", MMPAIndex.MMP_IGNORED)
			.executeUpdate();
		else
			throw new UserFriendlyException("Unsupported MMP query");
	}

	public void recalculatePairCounts() {
		logger.info("Recalculating pair counts...");
		Globals.session().createSQLQuery("update MMPTransformation set pairs_count=(select count(*) from MMPair where transformation_id=MMPTransformation.transformation_id)").executeUpdate();
		logger.info("Pair counts recalculated");
	}

	private long getRowsCount(String table) {
		return (Long) Globals.session().createSQLQuery("select count(*) c from " + table)
				.addScalar("c", LongType.INSTANCE)
				.uniqueResult();
	}

	@SuppressWarnings("unused")
	public void deleteInfrequentTransformations() {
		logger.info("Deleting infrequent transformations...");

		if(INFREQUENT_TRANSFORMATIONS <=0 )return;

		List<Long> ids = Globals.session().createQuery("select id from MMPTransformation where pairsCount < " + INFREQUENT_TRANSFORMATIONS).list();

		while (!ids.isEmpty())
		{
			List<Long> batchIds = ids.subList(0, Math.min(500, ids.size()));
			logger.info("Deleting " + batchIds.size() + " infrequent transformations out of " + ids.size());
			Globals.session().createSQLQuery("delete from MMPair where transformation_id in (:ids)").setParameterList("ids", batchIds).executeUpdate();
			Globals.session().createSQLQuery("delete from MMPTransformation where transformation_id in (:ids)").setParameterList("ids", batchIds).executeUpdate();
			batchIds.clear();
			Globals.restartAllTransactions(true);
		}

		MMPTransformation.clearCache();
		MMPAnnotationService.deleteInvalidAnnotations();
		MMPAnnotationService.countsCache.clear();

		logger.info("Infrequent transformations deleted");
	}

	public static void deleteInvalidPairs() {
		logger.info("Looking for invalid pairs associated with deleted transformations...");
		List<Long> pairsToDelete = Globals.session().createSQLQuery("select mmp_id from MMPair left join MMPTransformation t using (transformation_id) where t.transformation_id is null").addScalar("mmp_id", LongType.INSTANCE).list();
		logger.info("Deleting "+pairsToDelete.size()+" invalid pairs");
		while (!pairsToDelete.isEmpty()) {
			List<Long> batchIds = pairsToDelete.subList(0, Math.min(1000, pairsToDelete.size()));
			logger.info("Deleting " + batchIds.size() + " pairs");
			Globals.session().createSQLQuery("delete from MMPair where mmp_id in (:ids)").setParameterList("ids", batchIds).executeUpdate();
			batchIds.clear();
			Globals.restartAllTransactions(true);
		}
	}

	/**
	 * Removal of duplicated pairs coming to several transformations at once.
	 * This method using plan java code which may not work for mysql > 5.6
	 */
	public void deleteDublicatedPairs() {
		logger.info("Starting deletion of duplicated pairs by adding unique key...");

		String exception = "";

		try 
		{
			Class.forName(Globals.ochemConf.getProperty("hibernate.connection.driver_class")).newInstance();
			Connection conn = DriverManager.getConnection(Globals.ochemConf.getProperty("hibernate.connection.url"), Globals.ochemConf.getProperty("hibernate.connection.username"), Globals.ochemConf.getProperty("hibernate.connection.password"));
			Statement stat = conn.createStatement();
			stat.executeUpdate("alter ignore table MMPair add unique key un (mol1,mol2,transformation_id)");
			conn.close();

		} 
		catch (Exception e) 
		{
			exception = e.getMessage();
		} 

		try 
		{
			Class.forName(Globals.ochemConf.getProperty("hibernate.connection.driver_class")).newInstance();
			Connection conn = DriverManager.getConnection(Globals.ochemConf.getProperty("hibernate.connection.url"), Globals.ochemConf.getProperty("hibernate.connection.username"), Globals.ochemConf.getProperty("hibernate.connection.password"));
			Statement stat = conn.createStatement();
			stat.executeUpdate("alter table MMPair drop key un");
			conn.close();

		} 
		catch (Exception e) 
		{
			exception += e.getMessage();
		} 

		if(exception.length()>0)
			throw new UserFriendlyException("There were problems to delete duplicates because of "+exception);

		logger.info("Finished deletion of duplicated pairs by adding unique key...");
		Globals.restartAllTransactions(true);
	}

	/**
	 * Get the number of pending (init+assigned+ready) fragmentation tasks.
	 */
	private int getPendingTasksCount() throws MalformedURLException, IOException, ClassNotFoundException {
		int count = 0;
		List<Object[]> rows = new CalculationClient().getTasksSummary();
		for (Object[] row : rows)
		{
			String status = (String)row[1];

			if (QSPRConstants.MMPFrag.equals(row[0]))
				if (Task.isAliveStatus(status) || Task.READY.equals(row[1]) || Task.FETCHED.equals(row[1]))
					count += (Integer) row[2];
		}

		return count;
	}

	/**
	 * A cached value for forbidden transformation IDs
	 */
	private static Set<Long> forbiddenTransformationsIds = null;

	/**
	 * Retrieve the IDs of forbidden (ignored) transformations (e.g., annoying carbon chains).
	 * Pairs for such transformations should be ignored not to overload the database.
	 */
	public Set<Long> getForbiddenTransformationsIds() {

		if (forbiddenTransformationsIds != null)
			return forbiddenTransformationsIds;

		List<Long> ids = Globals.session().createSQLQuery("select transformation_id tid from MMPTransformation" + 
				"   inner join MMPFragment f1 on (fragment1_id=f1.fragment_id)" + 
				"   inner join MMPFragment f2 on (fragment2_id=f2.fragment_id)" + 
				"   where f1.carbon_chain and f2.carbon_chain and (f1.size > 1) and (f2.size > 1)")
				.addScalar("tid", LongType.INSTANCE)
				.list();

		Set<Long> uniqueIDs = new HashSet<Long>();
		uniqueIDs.addAll(ids);

		logger.info("Detected " + uniqueIDs.size() + " forbidden transformations");

		return forbiddenTransformationsIds = uniqueIDs;
	}

	private MMPIndexingService() {

	}

	private static MMPIndexingService instance = null;

	public static MMPIndexingService getInstance() {
		if (instance == null)
			instance = new MMPIndexingService();
		return instance;
	}

	public static void main(String[] args)
	{
		new WrapperThread()
		{

			@Override
			public void wrapped() throws Exception
			{
				//MMPIndexingService.getInstance().submitIndexingTasks();
				//System.out.println(MMPIndex.entryExists("aaa", 0L, 0));
				//MMPIndexingService.getInstance().fetchIndexingTasks();
				//MMPIndexingService.getInstance().reSubmitIndexingTasks(131627);
				//MMPAnnotationService.deleteInvalidAnnotations();
				MMPIndexingService.getInstance().rebuildMMPs();
				//MMPIndexingService.getInstance().deleteDublicatedPairs();
				//MMPIndexingService.getInstance().recalculatePairCounts();
				//MMPIndexingService.getInstance().deleteInfrequentTransformations(5);
				//MMPIndexingService.getInstance().deleteInvalidPairs();
			}
		}.run();
	}
}

