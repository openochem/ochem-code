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

package qspr.metaserver.cs;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.dao.MetalBondParserSdf;
import qspr.dao.Various;
import qspr.metaserver.ServerPool;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.cs.util.DescriptorsCacheTask;
import qspr.metaserver.cs.util.WebServiceCalculationTask;
import qspr.metaserver.util.MixtureAttachment;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.CalculationTask;
import qspr.workflow.utils.CalculationTaskSet;
import qspr.workflow.utils.QSPRConstants;
import qspr.metaserver.configurations.DescriptorType;

import com.eadmet.exceptions.CriticalException;

/*
 * Abstract descriptor calculator
 * Forwards the task to end-node servers (Dragon, EState, Fragmentor, workflow, ...)
 * 
 * @author: everyone has touched it
 */

public class DescriptorsServer extends WorkflowNodeServer
{
	private static transient final Logger logger = LogManager.getLogger(DescriptorsServer.class);

	private static final String CARBON = "\n\n\n  5  4  0  0  0  0  0  0  0  0999 V2000\n    1.0172   -0.0605   -0.0777 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.1094   -0.0605   -0.0777 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6532    0.6879    0.6296 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6532   -1.0473    0.2168 H   0  0  0  0  0  0  0  0  0  0  0  0\n    0.6532    0.1778   -1.0795 H   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  1  3  1  0  0  0  0\n  1  4  1  0  0  0  0\n  1  5  1  0  0  0  0\nM  END\n";

	/**
	 * Skip the cache completely, disregarding the task configuration
	 */
	public boolean isStandAlone;
	public boolean forceUpdateCache = false;

	public DescriptorsServer()
	{
		supportedTaskType = "Descriptors";
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);

		// For research purposes / Midnighter
		if (ServerPool.getInstance().getFreeServer(DescriptorsConfiguration.RANDOM) == null)
			ServerPool.getInstance().servers.add(new RandomGenerationServer());
	}

	public void setParam(String name, String value)
	{
		super.setParam(name, value);
		name = name.toUpperCase();
	}

	@Override
	public WorkflowNodeData calculate(WorkflowNodeData wndTask, Serializable receivedConfiguration) throws Exception
	{
		if (receivedConfiguration == null || !(receivedConfiguration instanceof DescriptorsConfiguration))
			throw new CriticalException("Invalid configuration passed, should be instance of DescriptorsConfiguration, got " + receivedConfiguration);
		DescriptorsConfiguration configuration = (DescriptorsConfiguration) receivedConfiguration;

		DataTable dtResult =  new DataTable(true);

		DataTable data = wndTask.ports.get(0);

		int initialData = data.getRowsSize(); // the size before expansion due to mixtures
		if (configuration.mixtures != null) {  // only if mixture processing is requested, we add single components at the end to get descriptors for them
			Set <String> existingMixtures = new HashSet<String>();
			for(int i = 0; i< initialData;i++)
				existingMixtures.add((String)data.getRow(i).getAttachment(QSPRConstants.INCHIKEYS));  // whatever we have

			for(int i = 0; i< initialData;i++) 
			{
				String componenst[] = Various.molecule.splitOrderedByChargeAndSize((String)data.getValue(i, 0));
				for(String s:componenst) {
					String inchi = Various.molecule.getInChIKeyNoStero(s);
					if(!existingMixtures.contains(inchi)){
						data.addRow();
						data.setValue(0, s); // storing sdf for calculation of descriptors
						data.getCurrentRow().addAttachment(QSPRConstants.INCHIKEYS, inchi);
						existingMixtures.add(inchi);
					}
				}
			}
		}


		if (configuration.types.size() == 0){
			out.println("No descriptor types specified!"); // can be OK for NoDescriptors
			for (int i = 0; i < data.getRowsSize(); i++)  // adding of attachments is done on the upper level and not required to be duplicated here
				dtResult.addRow().addAttachment(QSPRConstants.SDF_COLUMN, data.getValue(i, 0)); //adding 3D information for DIMENET
		}else{
			forceUpdateCache = (configuration.forceUpdateDescriptorCache != null && configuration.forceUpdateDescriptorCache);

			////forceUpdateCache  = true;

			List<DataTable> descTasks = calculateDescriptors(wndTask, configuration, true);
			if(descTasks == null)return new WorkflowNodeData(new DataTable(true));

			descTasks = selectDescriptors(descTasks, configuration,data);

			// Now create a combined answer
			dtResult = descTasks.get(0);
			dtResult.normalize();

			dtResult.id = "descriptors";

			out.println("Added descriptors (task 0) from " + configuration.types.get(0).type + ", " + dtResult.getColumnsSize() + " descriptors");

			boolean addPackageName = checkDuplicatedDescriptors(descTasks);

			if (addPackageName)
				out.println("Package names will be used to avoid duplicated descriptors.");

			if(configuration.allowsMissedValues())dtResult.allowNonStored();

			for (int i = 1; i < configuration.types.size(); i++)
			{
				DataTable dtNext = descTasks.get(i);
				dtNext.normalize(); 

				if(configuration.allowsMissedValues())dtNext.allowNonStored();

				String typeName = configuration.types.get(i).type;

				String name = addPackageName ? ":(" + typeName + ")" : "";
				int columns = typeName.equals(QSPRConstants.Workflow) ? 1 : dtNext.getColumnsSize();
				setStatus("Combining results from " + typeName + ", " + columns + " descriptors prefix: " + name);

				dtResult.mergeColumnsWith(dtNext, name); // adding suffix to merge columns, if required
			}

			if(configuration.allowsMissedValues())dtResult.markNoDescriptorsAsErrors();

		}

		if (configuration.mixtures != null)
			dtResult = MixtureDescriptorsProcessor.getInstance(configuration.mixtures).processMixtureDescriptors(data, dtResult); // will also get MixtureFractions from there

		// not done in browser to speed up
		for (int i = 0; i < dtResult.getRowsSize(); i++){ // store SMILES; required for NoDescriptors methods as well as for printing molecules identifiers
			AbstractDataRow row = dtResult.getRow(i);
			String smiles;
			try {
				String sdf = (String)data.getSDF(i);
				if(configuration.mixtures == null) // will be required later for correct validation
					row.addAttachment(QSPRConstants.MIXTURE_ATTACHMENT, new MixtureAttachment(sdf));

				if(row.getAttachment(QSPRConstants.SMILES_ATTACHMENT)==null) {// for mixtures already done
					smiles = Various.molecule.convertToFormatFixMetal(sdf, QSPRConstants.SMILESH);
					row.addAttachment(QSPRConstants.SMILES_ATTACHMENT, smiles); // Full kekulized SMILES without H atoms
				}

			}catch(Throwable e) { // do not add if there is an error
				System.out.println(e.getMessage());
				row.addAttachment(QSPRConstants.SMILES_ATTACHMENT, QSPRConstants.ERROR); // otherwise just report an error
				if(configuration.types.size() == 0)row.setError("Conversion to SMILES failed.");
				continue;
			}
		}

		dtResult = dtResult.getSlice(0, initialData);

		return new WorkflowNodeData(dtResult);
	}

	private List<DataTable> selectDescriptors(List<DataTable> descTasks, DescriptorsConfiguration configuration, DataTable original) throws Exception {
		if(configuration.supportSplicing()){
			DescriptorsConfiguration newConf = new DescriptorsConfiguration(configuration);
			newConf.types = new ArrayList<DescriptorType>(); 


			for(int i=0;i<configuration.types.size() ;i++)if(configuration.types.get(i).supportSplicing())
			{
				DescriptorType descType = configuration.types.get(i);
				descType.writeToCache = false; // we do not need to store 
				descType.skipCache = true; // we skip cache to have "fresh" calculations (in case if new descriptors appear, server was modified)
				newConf.types.add(descType);
			}

			DataTable data = new DataTable();
			data.addRow();
			data.setValue(QSPRConstants.SDF_COLUMN, CARBON);
			if(original.columnAttachments.get(QSPRConstants.SDF_COLUMN) != null)
				data.columnAttachments.put(QSPRConstants.SDF_COLUMN, original.columnAttachments.get(QSPRConstants.SDF_COLUMN));

			WorkflowNodeData wndTask = new WorkflowNodeData();
			wndTask.ports.add(data);

			List<DataTable> tab = calculateDescriptors(wndTask,newConf, false); // just for one molecule

			int j = 0;
			for(int i=0;i<configuration.types.size() ;i++)if(configuration.types.get(i).supportSplicing())
			{
				descTasks.get(i).keepByList(tab.get(j).getColumns());
				j++;
			}

		}
		return descTasks;
	}

	private List<DataTable> calculateDescriptors(WorkflowNodeData wndTask, DescriptorsConfiguration configuration, boolean calculateAllDescriptors) throws Exception {
		List<CalculationTask> descTasks = new ArrayList<CalculationTask>();

		if (wndTask.ports.get(0).getRowsSize() == 0){
			out.println("No molecules to calculate descriptors. Returning an empty datatable");
			return null;
		}

		for (DescriptorType descType : configuration.types){
			if(isStandAlone)descType.skipCache();
			logger.info("Processing " + descType.type);

			if(calculateAllDescriptors && descType.supportSplicing())
				descType = descType.configurationWithAllOnDescriptors();

			CalculationTask task = descType.type.startsWith(DescriptorsConfiguration.ApplyModel)
					?new WebServiceCalculationTask(descType, wndTask, this, forceUpdateCache)
							:new DescriptorsCacheTask(descType, wndTask, this, forceUpdateCache);
			task.post();
			descTasks.add(task);
		}

		CalculationTaskSet taskSet = new CalculationTaskSet(this);
		taskSet.tasks = descTasks;

		setStatus("Tasks are sent for calculations");
		taskSet.calculate(false);

		List<DataTable> results = new ArrayList<DataTable>();

		for(CalculationTask task : descTasks){

			if(task.getWndOutput() == null)throw new CriticalException("Task getWndOutput() is null in calculateDescriptors");
			if(task.getWndOutput().ports == null)throw new CriticalException("Task getWndOutput().ports is null in calculateDescriptors");
			if(task.getWndOutput().ports.get(0) == null)throw new CriticalException("Task getWndOutput().ports.get(0) is null in calculateDescriptors");

			results.add(task.getWndOutput().ports.get(0));
		}


		return results;
	}

	private boolean checkDuplicatedDescriptors(List<DataTable> descTasks)
	{
		HashSet<String> allDescriptors = new HashSet<String>();

		for (DataTable newDescriptors: descTasks)
			for (String descriptor : newDescriptors.getColumns())
			{
				if (allDescriptors.contains(descriptor))
					return true; // This column name was already there -- duplicates!!
				allDescriptors.add(descriptor);
			}

		return false;
	}

}
