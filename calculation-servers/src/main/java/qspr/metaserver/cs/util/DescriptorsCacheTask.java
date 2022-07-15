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

package qspr.metaserver.cs.util;

import java.io.IOException;
import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.interfaces.ProvidedConditions;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorType;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration.MixturesProcessing;
import qspr.metaserver.configurations.DescriptorsApplyModelConfiguration;
import qspr.metaserver.cs.MixtureDescriptorsProcessor;
import qspr.metaserver.cs.WorkflowNodeServer;
import qspr.metaserver.util.MixtureAttachment;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.utils.SDFProcessor;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.CalculationTask;

import com.eadmet.descriptorcache.DescriptorConfigEntry;
import com.eadmet.descriptorcache.DescriptorsCache;

public class DescriptorsCacheTask extends CalculationTask
{
	DescriptorType deskType;
	DataTable dsResult;  // store already retrieved descriptors from the cache

	private static transient final Logger logger = LogManager.getLogger(DescriptorsCacheTask.class);
	private static final String WAITING_FOR_RESULTS = "wait";

	public DescriptorsCacheTask(DescriptorType deskType, WorkflowNodeData wndTask, WorkflowNodeServer theParent, Boolean forceUpdateCache) throws Exception{

		parent = theParent;

		String user = parent == null ? null : parent.currentTask.get().getUser();
		//TODO deprecated code required for compatibility; Awkward!
		if(deskType.type.startsWith(DescriptorsConfiguration.ApplyModel)){  
			deskType.type = deskType.configuration.getDefaultTypeName();
			user = null;
		}

		taskName = deskType.type;
		configuration = deskType.configuration;
		this.deskType = deskType;

		logger.info("Preparing " + deskType.type);

		wndInput = wndTask.getDeeperCopy();

		if(configuration instanceof DescriptorsApplyModelConfiguration && wndInput.ports.size() > 1) { // we do not store descriptors with conditions (models) at the moment
			deskType.skipCache = true;
			deskType.writeToCache = false;
		}

		// Check for descriptors in storage
		if (!deskType.isSkipCache())
		{

			DescriptorConfigEntry dcEntry = new DescriptorConfigEntry(deskType);
			dcEntry.setUser(user); // required to possibly get specific descriptors stored by this user

			DescriptorsCache cache = new DescriptorsCache(deskType.markUncachedAsErrors == null ? false : deskType.markUncachedAsErrors, forceUpdateCache);

			// Check for descriptors in storage
			String condition = ProvidedConditions.isReplaceableCondition(deskType.type);
			if(deskType.type != null &&  condition != null && wndInput.ports.size()>1) { // only for stored solvent parameters
				DataTable mols = wndInput.ports.get(0).getDeepCopy();  // to eliminate problems with adding attachments
				DataTable cond = wndInput.ports.get(1);
				int initialSize = mols.getRowsSize();
				Map<String,Integer> ms = new HashMap<String,Integer>();
				for(int i = 0;i<initialSize;i++) {
					AbstractDataRow row = mols.getRow(i);
					MixtureAttachment ma = (MixtureAttachment)cond.getRow(i).getAttachment(QSPRConstants.SOLVENT_ATTACHMENT);
					if(ma != null) {
						if(row.getAttachment(QSPRConstants.MIXTURE_ATTACHMENT) != null)
						{
							row.setError("Currently solvent and mixtures cannot be used simultaneously");
							continue;
						}
						row.addAttachment(QSPRConstants.MIXTURE_ATTACHMENT, ma);
						int n=0;
						for(Map.Entry<String, Double> s:ma.fractions.entrySet()){
							String mol = ma.smiles.get(n);
							if(!ms.containsKey(mol)) {
								AbstractDataRow r = mols.addRow();
								ms.put(mol,mols.getRowsSize());
								r.setValue(0, "C");
								r.addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, mol);
								r.addAttachment(QSPRConstants.INCHIKEYS,s.getKey());
								//System.out.println("added: " + s.getKey() + " " + mols.currentRow);
							}
						}
					}else
						row.setError("No information about \"" + condition + "\" was provided");
				}

				//System.out.println("size mols " + mols.getRowsSize());

				WorkflowNodeData wnd = new WorkflowNodeData();
				wnd.setPort(0, mols);

				String originalType = null,userStored = user;

				if(dcEntry.type.contains(":")) {
					originalType = dcEntry.type;
					String parts[] = dcEntry.type.split(":");
					dcEntry.type = parts[1];
					dcEntry.setUser(user = parts[0]);
					dcEntry.updateMD5();
				}

				DataTable result = cache.getDescriptorsFromCache(dcEntry, wnd);

				if(originalType != null) {
					user = userStored;
					dcEntry.setUser(userStored); // restore initial entry
					dcEntry.type = originalType;
					dcEntry.updateMD5();
				}

				MixtureDescriptorsProcessor fraction = MixtureDescriptorsProcessor.getInstance(MixturesProcessing.SUMONLY);
				result = fraction.processMixtureDescriptors(mols, result);

				result = result.getSlice(0, initialSize);

				for(int i =0; i<result.getRowsSize();i++) {
					result.setValue(i, condition + QSPRConstants.USED, 1.); // to prevent using condition on the selection step
					if(result.getRow(i).getAttachment(QSPRConstants.MIXTURE_ATTACHMENT) != null)
						result.getRow(i).attachments.remove(QSPRConstants.MIXTURE_ATTACHMENT); // to remove these MIXTURES....

					if(result.getRow(i).isError())
						result.getRow(i).setError("No stored descriptors were found for (one of) " + condition);

					for(int j =0; j<result.getColumnsSize();j++)
						if(((Double)result.getValue(i, j)).isNaN()) {
							result.getRow(i).setError("Some stored values were NaN at least for colum " + result.getColumn(j));
							break;
						}
				}

				if(result.getRowsNoErrorsSize() == 0)
					throw new Exception("No stored descriptors were found for all records: check descriptor type: " + deskType.type);
				else
					System.out.println("Found actual solvents for " + result.getRowsNoErrorsSize());


				dsResult = addNewlyCalculated(result);
				wndOutput = new WorkflowNodeData(dsResult);  // work is done!
			}
			else {
				while(wndInput.ports.size() > 1) // These descriptors do not need any conditions; remove the conditions
					wndInput.ports.remove(1);

				dsResult = cache.getDescriptorsFromCache(dcEntry, wndInput);

				eliminateDuplicates(wndInput, dsResult, dcEntry);

				if(dsResult.getRowsSize() != wndInput.getRows()) throw new Exception("not equal number of rows");

				DataTable mols = wndInput.ports.get(0);
				DataTable dtUncached = mols.getEmptyCopy();

				for(int i=0;i<dsResult.getRowsSize();i++)
					if(DescriptorsCache.SENT_FOR_CALCULATIONS.equals(dsResult.getRow(i).status))
						dtUncached.addRow(mols.getRow(i));
				wndInput.setPort(0, dtUncached);

				// Anything left for real calculation?
				if (dtUncached.getRowsSize() > 0){
					if(configuration instanceof DescriptorsAbstractConfiguration){
						DescriptorsAbstractConfiguration conf = (DescriptorsAbstractConfiguration)configuration;
						deskType.writeToCache = conf.isCachable(); // this particular configuration is not cachable
						if(conf.isLongCalculation() && dtUncached.getRowsSize() < 1000 && dtUncached.getRowsSize() < dsResult.getRowsSize()/10)
							conf.repostSize = (conf.repostSize == null || conf.repostSize/10 == 0)? 1 : conf.repostSize/10; // typical scenario - calculation of few molecules that are failing and blocking everything
					}
				}
				else
					wndOutput = new WorkflowNodeData(dsResult);  // work is done!
			}
		}
	}

	private void eliminateDuplicates(WorkflowNodeData wndInput, DataTable dsResult, DescriptorConfigEntry dcEntry) {

		DataTable mol = wndInput.ports.get(0);
		Map<String,Integer> sentData = new HashMap<String, Integer>(); // are used to avoid duplicated structures

		for(int n=0;n<mol.getRowsSize();n++){

			if(!DescriptorsCache.SENT_FOR_CALCULATIONS.equals(dsResult.getRow(n).status))continue;

			String sdf = (String) mol.getValue(n,0);

			String moleculeMD5 = sdf != null ? SDFProcessor.getMD5SDF(sdf,dcEntry.twoD) : "null";

			if(sentData.containsKey(moleculeMD5)){
				dsResult.getRow(n).setStatus(WAITING_FOR_RESULTS + " " + sentData.get(moleculeMD5)); // there is no need to send this molecule again
				continue;
			}

			sentData.put(moleculeMD5, n); // the first entry
		}		
	}

	@Override
	public void post() throws ClassNotFoundException, IOException, InterruptedException{
		logger.info("Starting " + deskType.type);
		if(wndOutput == null)super.post();
	}

	@Override
	public void onReady() throws Exception
	{
		DataTable newValues = wndOutput.ports.get(0);  // something coming from calculations

		if (deskType.writeToCache != null && deskType.writeToCache){
			DescriptorsCache cache = new DescriptorsCache();
			DescriptorConfigEntry dcEntry = new DescriptorConfigEntry(deskType);
			cache.saveNewValues(dcEntry, wndInput, newValues);
		}

		wndOutput.ports.set(0, addNewlyCalculated(newValues));
	}

	/**
	 * Fill in rows with newly calculated descriptors
	 * @param newDescriptors
	 * @throws IOException 
	 */

	private DataTable addNewlyCalculated(DataTable newDescriptors) throws IOException {

		if(dsResult == null) return newDescriptors;

		DataTable dtComposedResult = new DataTable(true);

		boolean newIsLarger = newDescriptors.getRowsSize() > dsResult.getRowsSize()/2;

		dtComposedResult.addColumns(newIsLarger ? newDescriptors.getColumns() : dsResult.getColumns());  // bigger set will be just copied without creating new rows

		for (String column : newIsLarger ? dsResult.getColumns() : newDescriptors.getColumns())
			if(!dtComposedResult.containsColumn(column))
				dtComposedResult.addColumn(column);

		dtComposedResult.addRowsFrom(dsResult);
		dtComposedResult.addRowsFrom(newDescriptors); // we first add newDescriptors at the end of the list

		int n = dsResult.getRowsSize();
		for(int i=0; i<dsResult.getRowsSize(); i++){
			String status = dsResult.getRow(i).status;
			if(status == null)status="";

			if(status.equals(DescriptorsCache.SENT_FOR_CALCULATIONS)) // was sent for calculation
				dtComposedResult.setRow(i,dtComposedResult.getRow(n++)); // and substitute those sent for calculations with new values
			else
				if(status.startsWith(WAITING_FOR_RESULTS)){
					Integer row = Integer.valueOf(status.substring(WAITING_FOR_RESULTS.length()).trim());
					Map<String, Serializable> attachments = dtComposedResult.getRow(i).attachments;
					dtComposedResult.setRow(i,dtComposedResult.getRow(row).getDeeperCopy());
					dtComposedResult.getRow(i).attachments = attachments;
				}		
		}

		dtComposedResult.setSize(dsResult.getRowsSize());

		return dtComposedResult;
	}


}
