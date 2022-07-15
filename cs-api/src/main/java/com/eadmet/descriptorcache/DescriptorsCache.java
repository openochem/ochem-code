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

package com.eadmet.descriptorcache;

import java.io.Serializable;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.utils.SDFProcessor;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

/**
 * A service for managing the descriptors storage
 * @author midnighter
 *
 */
public class DescriptorsCache
{

	private static final Logger logger = LogManager.getLogger(DescriptorsCache.class);

	final static public String SENT_FOR_CALCULATIONS ="sent";
	final static public int BATCH = 10000; // number of descriptors downloaded in one batch; implemented to decrease memory for memory expensive cache records

	/**
	 * When set to "true", the cache will mark even properly cached descriptors for recalculation
	 */

	boolean forceUpdateCache = false;

	/**
	 * We need to distinguish stored and non-stored descriptors for their proper merging
	 */
	boolean storedDescriptors = false;

	public DescriptorsCache (boolean storedDesciptors, boolean forceUpdateCache )
	{
		this.storedDescriptors = storedDesciptors;
		this.forceUpdateCache = forceUpdateCache;
	}

	public DescriptorsCache()
	{
	}

	public void saveCacheEntries(DescriptorConfigEntry dcEntry, List<CacheEntry> entries) throws Exception
	{
		DescriptorsRepository repository = DescriptorsRepositoryFactory.getReattemptingRepository();
		String[] md5 = new String[entries.size()];

		for (int i = 0; i < md5.length; i++)
		{
			md5[i] = entries.get(i).moleculeMD5;
			if (!entries.get(i).config.equals(dcEntry))
				throw new RuntimeException("Cannot mix configurations while saving cache entries");
		}

		// Query for duplicates
		List<CacheEntry> savedEntries = repository.getDescriptors(md5, dcEntry);

		// Set objectIDs for the duplicates so that they are overwritten
		for (int i = 0; i < savedEntries.size(); i++)
			if (savedEntries.get(i) != null)
				entries.get(i).objectID = savedEntries.get(i).objectID;

		repository.saveDescriptors(entries);
	}

	/**
	 * Retrieve descriptors from the cache
	 * If descriptors are absent, return empty rows 
	 * @param dtMolecules
	 * @param descriptorType
	 * @param configurationXML
	 * @return
	 * @throws Exception
	 */

	public DataTable getDescriptorsFromCache(DescriptorConfigEntry dcEntry, WorkflowNodeData wndInput) throws Exception{

		DataTable dtResults = new DataTable(true);
		Map<String, Integer> columns = new HashMap<String, Integer>(); // are used to speed up search in the list

		List<CacheEntry> cacheEntries = requestEntriesBatch(dcEntry, wndInput, BATCH);

		for (CacheEntry cacheEntry : cacheEntries){

			dtResults.addRow();

			if (cacheEntry.sendForCalculation())
				dtResults.getCurrentRow().setStatus(SENT_FOR_CALCULATIONS);
			else
			{
				if (cacheEntry.error != null)
					dtResults.getCurrentRow().setError("Calculation failed after " + cacheEntry.attempts + " attemts. Exemplary error: " + cacheEntry.error);
				else
				{
					// A valid cached result
					String[] names = cacheEntry.getNames();
					float[] values = cacheEntry.getValues();

					for (int i = 0; i < names.length; i++)
					{
						Integer index = columns.get(names[i]);
						if (index == null)
						{
							dtResults.addColumn(names[i]);
							columns.put(names[i], index = dtResults.getColumnsSize() - 1);
						}

						dtResults.setValue(index, values[i]);
					}
				}
			}
		}

		return dtResults;

	}

	private List<CacheEntry> requestEntriesBatch(DescriptorConfigEntry dcEntry, WorkflowNodeData wndInput, int batchSize) throws Exception{

		List<CacheEntry> list = new ArrayList<CacheEntry>();
		int molSize = wndInput.getRows();
		DataTable mols, cond = null;

		for(int mol=0; mol< molSize; mol += batchSize){

			logger.info("Processing " + mol + " out of " + molSize + " using batch " + batchSize);

			mols = wndInput.ports.get(0).getSlice(mol, mol + batchSize > molSize ? molSize : mol + batchSize );
			if(wndInput.ports.size() >1)
				cond =  wndInput.ports.get(1).getSlice(mol, mol + batchSize > molSize ? molSize : mol + batchSize );

			list.addAll(requestEntries(dcEntry, mols, cond));
		}

		return list;
	}


	/**
	 * Request cache entries for the molecules and given Entry
	 * Consider molecule MD5, MP2, and external IDs stored in the dataTable attachments
	 * @return 
	 * 
	 */

	private List<CacheEntry> requestEntries(DescriptorConfigEntry dcEntry, DataTable dtMolecules, DataTable conditions) throws Exception{
		List<CacheEntry> cacheEntriesSdf = requestEntriesSDF(dcEntry, dtMolecules, conditions, true);
		List<CacheEntry> cacheEntriesMD5 = requestEntriesSDF(dcEntry, dtMolecules, conditions, false);

		for(int i = 0; i< cacheEntriesSdf.size(); i++)
			if(cacheEntriesSdf.get(i).error != null && cacheEntriesMD5.get(i).error == null) // something was found ...
				cacheEntriesSdf.set(i,cacheEntriesMD5.get(i));
		return cacheEntriesSdf;
	}

	private List<CacheEntry> requestEntriesSDF(DescriptorConfigEntry dcEntry, DataTable dtMolecules, DataTable conditions, boolean useSDF) throws Exception
	{
		List<CacheEntry> cacheEntries;

		DescriptorsRepository repository = DescriptorsRepositoryFactory.getReattemptingRepository();
		String user = dcEntry.user; dcEntry.setUser(null);  //User information is not required at this point
		repository.saveConfig(dcEntry);

		String[] md5 = new String[dtMolecules.getRowsSize()];
		Integer[] mp2 = new Integer[dtMolecules.getRowsSize()];
		dtMolecules.reset();
		while (dtMolecules.nextRow())
		{
			String externalId = (String) dtMolecules.getCurrentRow().getAttachment(QSPRConstants.EXTERNAL_ID);
			
			String sdf = useSDF || externalId == null? (String) dtMolecules.getValue() : null;
			if (sdf != null && sdf.length() > 0)
			{
				md5[dtMolecules.currentRow] = SDFProcessor.getMD5SDF(sdf,dcEntry.twoD);
				Serializable attachment = dtMolecules.getCurrentRow().getAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM);
				mp2[dtMolecules.currentRow] = attachment != null ? Integer.valueOf("" + attachment) : -1;
				if(mp2[dtMolecules.currentRow] < 0)
					mp2[dtMolecules.currentRow] = null; // MP2 should be only positive ones. Otherwise they are temporal ones and should not be used! 
			}
			else
			{
				md5[dtMolecules.currentRow] = externalId;
				mp2[dtMolecules.currentRow] = null;
			}
		}

		// getting all public data for User null first
		cacheEntries = repository.getDescriptors(md5, dcEntry);

		// Override with published entries if any
		dcEntry.setUser(QSPRConstants.PUBLISHER);
		List<CacheEntry> privateEntries = repository.getDescriptors(md5, dcEntry); //Private by md5
		updateEntries(privateEntries, cacheEntries);

		dcEntry.setUser(user); // restore initial entry

		// Override with private entries if any
		if (user != null)
		{
			privateEntries = repository.getDescriptors(md5, dcEntry); //Private by md5
			updateEntries(privateEntries, cacheEntries);

			privateEntries = repository.getDescriptors(mp2, dcEntry); //Private by mp2
			updateEntries(privateEntries, cacheEntries);
		}

		int cntNew = 0;
		int cntRetryErrors = 0;
		for (int i = 0; i < cacheEntries.size(); i++)
		{
			if (cacheEntries.get(i) == null)
			{
				CacheEntry entry = new CacheEntry();
				entry.dateCreated = new Timestamp(Calendar.getInstance().getTimeInMillis());
				entry.setMD5(md5[i]);
				entry.mp2 = mp2[i];
				entry.config = dcEntry;
				if (storedDescriptors || (dtMolecules.getValue(i, 0) == null))
				{
					entry.doNotSend();
					entry.error = QSPRConstants.NOSTORED + (entry.objectID == null?"":(" for " + entry.objectID));
					cacheEntries.set(i, entry);
				}
				else
				{
					entry.isNew = true;
					cacheEntries.set(i, entry);
					if (cacheEntries.get(i).sendForCalculation())
						cntNew++;
				}
			}
			else
			{
				if (forceUpdateCache)
				{
					cacheEntries.get(i).doSend();
					cacheEntries.get(i).attempts = 0;
				}

				if (cacheEntries.get(i).sendForCalculation())
					cntRetryErrors++;
			}
		}

		logger.info(String.format("%d cached entries, %d will be sent for calculation", dtMolecules.getRowsSize() - cntNew, cntNew + cntRetryErrors));
		if (cntRetryErrors != 0)
			if (forceUpdateCache)
				logger.info("Force recalculation of " + cntRetryErrors + " cache entries will be performed");
			else
				logger.info("Calculation of " + cntRetryErrors + " error entries will be re-attempted");
		return cacheEntries;
	}


	void updateEntries(List<CacheEntry> privateEntries, List<CacheEntry> cacheEntries){
		for (int i = 0; i < privateEntries.size(); i++)
			if (privateEntries.get(i) != null && !privateEntries.get(i).sendForCalculation()) // always prefer to use global cache
				cacheEntries.set(i, privateEntries.get(i));
	}

	/**
	 *  Saves new calculated values for molecules
	 * @param dcEntry
	 * @param dtMolecules
	 * @param dtResults
	 * @throws Exception
	 */

	public synchronized void saveNewValues(DescriptorConfigEntry dcEntry, WorkflowNodeData wndInput, DataTable dtResults) throws Exception
	{
		if(wndInput.getRows() != dtResults.getRowsSize()) throw new Exception("Different sizes for wndInput.getRows() != dtResults.getRowsSize() " + 
				wndInput.getRows() + " != " + dtResults.getRowsSize());

		String[] names = dtResults.getColumns().toArray(new String[0]);
		float[] values = new float[names.length];

		List<CacheEntry> cacheEntries = requestEntriesBatch(dcEntry, wndInput, BATCH);
		List<CacheEntry> newEntries = new ArrayList<CacheEntry>();

		Set<String> alreadyHere = new HashSet<String>();

		for(int j=0; j< cacheEntries.size(); j++){
			CacheEntry entry = cacheEntries.get(j);
			if (!entry.sendForCalculation() && entry.error == null)continue; // this entry was already saved possibly by another server and it is not an error 

			if(alreadyHere.contains(entry.moleculeMD5)) continue;
			alreadyHere.add(entry.moleculeMD5);

			newEntries.add(entry);

			AbstractDataRow row = dtResults.getRow(j); // getting respective calculated result

			if (row.isError())
			{
				entry.error = row.detailedStatus;
				entry.attempts++;
			}
			else
			{
				entry.error = null;
				for (int i = 0; i < values.length; i++)
					if(names[i].contains(QSPRConstants.INDIVIDUAL_PREDICTIONS))
						values[i] = 0;
					else
						values[i] = Float.valueOf("" + row.getValue(i));
				entry.setNamesAndValues(names, values, dcEntry.keepAllValues);
			}
		}

		logger.info(String.format("Saving %d new cache entries", newEntries.size()));

		DescriptorsRepository repository = DescriptorsRepositoryFactory.getReattemptingRepository();
		if(newEntries.size()>0)
			repository.saveDescriptors(newEntries);
	}


	public static void main(String[] args) throws InterruptedException
	{
	}

	public static void deleteCachedEntries(DescriptorConfigEntry dcEntry,
			String user) throws Exception {
		DescriptorsRepository repository = DescriptorsRepositoryFactory.getReattemptingRepository();
		repository.saveConfig(dcEntry);
		repository.clearCache(dcEntry);		
	}
}
