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

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Calendar;
import java.util.List;

import javax.xml.bind.JAXBException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.eadmet.descriptorcache.DescriptorConfigEntry;
import com.eadmet.descriptorcache.DescriptorsCache;
import com.eadmet.exceptions.UserFriendlyException;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ConditionSet;
import qspr.entities.Model;
import qspr.metaserver.configurations.DescriptorsApplyModelConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

/**
 * The class stores cached values as well as experimental values for models
 * If model is published, values will be cached
 */


public class ModelApplierCache implements Serializable
{
	private static transient final Logger logger = LogManager.getLogger(ModelApplierCache.class);

	private static final long serialVersionUID = 1L;

	public static final String LOSTRESULT = "Cached result was lost";

	// Required because of incompatibility with ApplyModel results (different columns due to conversion of units)
	private static final String STORED_PREDICTIONS = "StoredPredictions";

	public transient Basket basket;
	public transient Model model;
	ConditionSet defaultConditions;

	// Introduces memory overheads to store records in memory
	private transient DataTable storedPredictions;

	BitSet cachedEntries = new BitSet(); // on indicate that the corresponding entry is available as a cached entry

	public ModelApplierCache(){
	}

	/**
	 *  Creates model 
	 * @param model
	 * @param basket
	 */

	public ModelApplierCache(Model model){
		this.model = model;
	}

	public boolean isEmpty(){
		return getCachedCount() == 0;
	}

	public ModelApplierCache(Model model, Basket basket, PredictionScenario type, ConditionSet defaultConditions){
		this.basket = basket;
		this.model = model;
		this.defaultConditions = defaultConditions;

		// 			ModelProcessor processor = ModelFactory.getProcessor(model);
		//		task.setConfiguration(processor.getApplierConfiguration(recalculated));
		//	WorkflowNodeData wn = (WorkflowNodeData) processor.getDataForApplier(basket, parent.defaultConditions);

		DataTable stored = getStoredPredictions(true);
		for (int i=0; i < basket.entries.size(); i++)
			if(!DescriptorsCache.SENT_FOR_CALCULATIONS.equals(stored.getRow(i).status) ) // was found!
				cachedEntries.set(i); // this entry is cached

		logger.info(getCachedCount() + " unique molecules are present in cache and won't be sent for calculation");
	}



	public Basket getNonCachedEntriesBasket(){
		Basket newb = new Basket();
		newb.entries = getNonCachedEntries(basket.entries);
		return newb;
	}

	private DataTable getStoredPredictions(boolean all){
		Basket stored = new Basket();
		DataTable res = null;

		if(storedPredictions != null) {
			res = new DataTable(true);
			res.columns = storedPredictions.columns;
		}

		for(int i = 0 ; i <basket.entries.size();i++)
			if(all || cachedEntries.get(i)){ // was cached and was not sent for calculations or all entries are required
				stored.entries.add(basket.entries.get(i));
				if(storedPredictions != null) res.addRow(storedPredictions.getRow(i));
			}

		if(storedPredictions != null) return res;

		fetchMoleculesForBasketEntries(stored);

		try{
			ModelProcessor processor = ModelFactory.getProcessor(model);
			//		task.setConfiguration(processor.getApplierConfiguration(recalculated));
			WorkflowNodeData wn = (WorkflowNodeData) processor.getDataForApplier(stored, defaultConditions);
			DescriptorsCache cache = new DescriptorsCache(false, false);
			DescriptorConfigEntry dcEntry = createConfig(model.publicId);
			DataTable tab = cache.getDescriptorsFromCache(dcEntry, wn); // get cached records
			DataTable mol = wn.ports.get(0);
			for(int i = 0; i< tab.getRowsSize(); i++){
				String sdf = (String) mol.getValue(i,0);
				if(sdf == null || sdf.length() == 0)
					tab.getRow(i).setStatus(DescriptorsCache.SENT_FOR_CALCULATIONS); // These molecules can be based on external IDs
				else {
					try {
						tab.getRow(i).addAttachment(QSPRConstants.SMILES_ATTACHMENT, Various.molecule.convertToCanonicalSMILES( (String)sdf));
					}catch(IOException e) {
						tab.getRow(i).addAttachment(QSPRConstants.SMILES_ATTACHMENT, QSPRConstants.ERROR_SMILES);
					}
					String inchies = Various.molecule.getInChiKey(sdf);
					if(inchies.contains("-"))
					tab.getRow(i).addAttachment(QSPRConstants.INCHIKEYS, 
							Various.molecule.getInChiKey(sdf)
							);
				}
			}

			tab.reset();
			while(tab.nextRow())
				if(!DescriptorsCache.SENT_FOR_CALCULATIONS.equals(tab.getCurrentRow().status))
					tab.getCurrentRow().addAttachment(QSPRConstants.CACHED, true);

			if(all)storedPredictions = tab;
			return tab;
		}catch(Exception e){

			DataTable mols = new DataTable();
			for(int i =0; i<stored.entries.size();i++ )
				mols.addRow().setError(e.getMessage());
			return mols;
		}
	}

	private void fetchMoleculesForBasketEntries(Basket basket){
		for(BasketEntry entry: basket.entries)
			entry.ep.molecule = Repository.molecule.getMolecule(entry.ep.molecule.id);
	}

	private static DescriptorConfigEntry createConfig(Long publicId) throws JAXBException{
		DescriptorsApplyModelConfiguration config = new DescriptorsApplyModelConfiguration(publicId);
		return new DescriptorConfigEntry(config, STORED_PREDICTIONS + "_" + publicId);
	}

	private List<BasketEntry> getNonCachedEntries(List<BasketEntry> originalList){
		List<BasketEntry> entries = new ArrayList<BasketEntry>(originalList);
		for(int i=entries.size()-1;i>=0;i--)
			if(cachedEntries.get(i))
				entries.remove(i);
		return entries;
	}

	public void cachePredictions(DataTable dtResults)
	{
		logger.info("Profiling: Starting caching..."); 

		long start = Calendar.getInstance().getTimeInMillis();

		final Basket newb = new Basket();
		newb.entries = 	getNonCachedEntries(basket.entries);
		newb.evict(true);
		final DataTable tab = dtResults.getCopy(); // new copy is required for caching; dtResults will be changed in above thread

		logger.info("Profiling: Fetching normally... " + newb.entries.size()); 

		Thread t = new Thread()
		{
			@Override
			public void run()
			{
				long start = Calendar.getInstance().getTimeInMillis();
				Globals.startAllTransactions();
				try{
					model = Model.initialiseModel(model);
					logger.info("Profiling: Entering asyncronous mode for " + model.publicId); 
					fetchMoleculesForBasketEntries(newb); // get molecule for each entry
					logger.info("Profiling: Fetching molecules in asyncronous mode took " + (Calendar.getInstance().getTimeInMillis() - start)+"ms");  start = Calendar.getInstance().getTimeInMillis();
					ModelProcessor processor = ModelFactory.getProcessor(model);
					WorkflowNodeData wn = (WorkflowNodeData) processor.getDataForApplier(newb, defaultConditions);
					saveCache(wn, tab);
				}catch(Exception e){
					e.printStackTrace();
					logger.error("Profiling:  " + e.getMessage() + " after " + (Calendar.getInstance().getTimeInMillis() - start)+"ms");  start = Calendar.getInstance().getTimeInMillis();
				}
				Globals.rollbackAllTransactions();
				logger.info("Profiling: Saving caches in asyncronous mode took " + (Calendar.getInstance().getTimeInMillis() - start)+"ms for model " + model.publicId);  start = Calendar.getInstance().getTimeInMillis();
			}
		};
		t.start();

		logger.info("Profiling: Finished " + (Calendar.getInstance().getTimeInMillis() - start)+"ms"); 

	}

	/**
	 * Saves new values to the cache
	 * @param dtPredictions
	 * @param mols
	 * @param type
	 * @throws Exception 
	 */

	private void saveCache(WorkflowNodeData wn, DataTable dtResults){

		if(dtResults.getRowsNoErrorsSize() == 0) return; // nothing to be cached

		int molSize = wn.getRows();

		if(molSize != dtResults.getRowsSize()) 
			throw new UserFriendlyException("Different numbers of rows in saveCache " + molSize +" != " + dtResults.getRowsSize());

		try{
			for(int i=0; i< dtResults.getRowsSize();i++)
				if(dtResults.getRow(i).isError()){
					for(DataTable tab : wn.ports)
						tab.deleteRow(i);
					dtResults.deleteRow(i);
				}

			DescriptorConfigEntry dcEntry = createConfig(model.publicId);
			DescriptorsCache cache = new DescriptorsCache(false, false);
			cache.saveNewValues(dcEntry, wn, dtResults);
			logger.info("" + molSize + " new predictions saved to the cache");
		}catch(Exception e){
			logger.error("saveCache failed failed: " + e.getMessage());
		}
	}


	/**
	 * Get all results -- either everything was cached or calculation task failed
	 * @param errorMsgForMissingRows
	 * @return
	 */

	public DataTable getCachedResult(String errorMsgForMissingRows){
		return getCachedResult(errorMsgForMissingRows, true);
	}

	private DataTable getCachedResult(String errorMsgForMissingRows, boolean allRecords){
		DataTable cachedResult = getStoredPredictions(allRecords); // if there are no new results, we will return all entries 	
		cachedResult.reset();
		while(cachedResult.nextRow())
			if(DescriptorsCache.SENT_FOR_CALCULATIONS.equals(cachedResult.getCurrentRow().status))
				cachedResult.getCurrentRow().setError(errorMsgForMissingRows); // marking non processed records
		return cachedResult;
	}

	public DataTable mergeCachedResult(DataTable newResults){

		if(cachedEntries.cardinality() == 0)return newResults; //nothing has been yet cached!

		DataTable cachedResult = getCachedResult(LOSTRESULT,false);
		if(newResults == null) return cachedResult;

		if(newResults.getRowsNoErrorsSize() == 0) // all new predictions failed
			newResults.columns = cachedResult.columns; 

		if(cachedResult.getRowsNoErrorsSize() == 0) // all cached results are errors
			cachedResult.columns = newResults.columns; 

		if(newResults.compareColumns(cachedResult) != null && newResults.getRowsNoErrorsSize() > 0)
			for(int i = 0; i < cachedResult.getRowsSize() ; i++)
				cachedResult.getRow(i).setError("Cached and calculate columns are different: clean cache results! " + newResults.compareColumns(cachedResult));

		int n = 0;
		for(int i=0; i < basket.entries.size(); i++)
			if(cachedEntries.get(i))
				newResults.insertRow(i, cachedResult.getRow(n++)); // adding cached results

		return newResults;
	}

	public int getCachedCount(){
		return cachedEntries == null ? 0 : cachedEntries.cardinality();
	}

	public static void clearCachedPredictions(Long publicId) {
		try{
			DescriptorConfigEntry entry = createConfig(publicId);
			DescriptorsCache.deleteCachedEntries(entry,null); // Deleting stored descriptors
			DescriptorsApplyModelConfiguration conf = new DescriptorsApplyModelConfiguration(publicId);
			entry.type = conf.getDefaultTypeName(); // Deleting ApplyModel descriptors
			DescriptorsCache.deleteCachedEntries(entry,null);
		}catch(Exception e){
			System.out.println(e.getMessage());
		}
	}

}

