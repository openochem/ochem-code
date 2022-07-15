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

package qspr.metaserver.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import qspr.metaserver.configurations.ValidationConfiguration;
import qspr.metaserver.cs.ValidationAbstractServer;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.CalculationTask;
import qspr.workflow.utils.QSPRConstants;

public class MachineLearningTask extends CalculationTask
{
	public Integer trainingSetSize = 0;
	public List<Integer> validationRowNums = new ArrayList<Integer>();
	public Map<Integer, Integer> trainingRowRepetitions = new HashMap<Integer, Integer>();

	/**
	 * Indicates whether this model has to be kept for correct work of the Validation Protocols
	 * It has nothing to do whether we want to save resulting models!
	 */
	public boolean keepModel = true;

	/**
	 * 
	 * @param server
	 * @param wndOriginal
	 * @param validationConfiguration
	 * @param net  if 
	 * type == configuration.ensembleSize, all molecules will be used
	 * @throws IOException 
	 */
	public MachineLearningTask(ValidationAbstractServer server, WorkflowNodeData wndOriginal, ValidationConfiguration validationConfiguration, int net)
			throws IOException
	{
		setParent(server);
		client.out = server.out;
		client.setSid("Validation");

		setId("analysis " + net);

		// Copy data table structures from original task
		wndInput = wndOriginal.getEmptyCopy();
		configuration = validationConfiguration.taskConfiguration;
		taskName = validationConfiguration.taskName;
	}

	@Override
	public void onReady()
	{
		if (!keepModel)
		{
			if (wndOutput != null && wndOutput.ports != null && wndOutput.ports.size() > 1)
			{
				wndOutput.ports.set(1, new DataTable());
				client.out.println("Removing model configuration (its not needed, just to save memory)");
			}
		}
		System.gc();
	}

	/**
	 * Create model that is cope of the training set
	 * @param wndOriginal
	 * @param splitter
	 * @param addValidationSet
	 */

	public void replicateTrainingSet(WorkflowNodeData wndOriginal)
	{
		if(wndOriginal.ports.size() == 1) wndOriginal.ports.add(new DataTable()); // FIX to have a possibility of a normal apply of models
		trainingSetSize = wndOriginal.ports.get(1).getRowsSize();

		Set<Integer>  inTraining = new HashSet<Integer>();

		for (int record = 0; record < wndOriginal.ports.get(0).getRowsSize(); record++)
		{ // to preserve the same order as in the provided set
			AbstractDataRow r = wndOriginal.ports.get(0).getRow(record);
			Integer rid = (Integer)r.getAttachment(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM);
			if (record < trainingSetSize) {
				addTrainingRecord(record, 1, wndOriginal);
				inTraining.add(rid);
			}
			else {
				boolean exclude = r.getAttachment(QSPRConstants.EXCLUDED) != null;
				if(!exclude || !inTraining.contains(rid)) // not excluded or not part of the training set
					addValidationRecord(record,r);
			}
		}
	}

	/**
	 * 
	 * @param trainingSetMoleculeHashesRepetitons
	 * @param wndOriginal
	 * @param splitter
	 * @param fullSetTask  -- indicate to add to the task also external validation set -- it is not required for individual cross-validation loops
	 */
	public void createTrainingAndValidationSets(Map<Integer, Integer> trainingSetMoleculeHashesRepetitons, WorkflowNodeData wndOriginal,
			TrainingSetSplitter splitter, boolean fullSetTask)
	{
		Set<Integer> trainingSetRecords = new TreeSet<Integer>();

		Set<Integer>  inTraining = new HashSet<Integer>(), inValidation = new HashSet<Integer>();

		for (Integer molecule: trainingSetMoleculeHashesRepetitons.keySet())
			for (Integer record: splitter.getRecords(molecule))
			{
				if(trainingSetRecords.contains(record))continue; // already added
				addTrainingRecord(record, trainingSetMoleculeHashesRepetitons.get(molecule), wndOriginal);
				trainingSetRecords.add(record);
				inTraining.add((Integer)wndOriginal.ports.get(0).getRow(record).getAttachment(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM));
			}

		trainingSetSize = wndInput.ports.get(1).getRowsSize();

		int originalTrainingSetSize = wndOriginal.ports.get(1).getRowsSize();

		// add all other non-selected records as the validation set
		for (int record = 0; record < wndOriginal.ports.get(0).getRowsSize(); record++) {
			AbstractDataRow r = wndOriginal.ports.get(0).getRow(record);
			Integer rid = (Integer)r.getAttachment(QSPRConstants.MOLECULE_ID_NO_STEREOCHEM);

			if (record < originalTrainingSetSize )
			{
				if(!trainingSetRecords.contains(record)) {
					trainingSetRecords.add(record); // to avoid adding duplicates from that set
					addValidationRecord(record,r);
					inValidation.add(rid);
				}
				continue;
			}

			boolean exclude = r.getAttachment(QSPRConstants.EXCLUDED) != null;

			if(fullSetTask) { // we add only records that are not EXCLUDED or those that are not part of the full training set
				if(!exclude || !inTraining.contains(rid)) // not excluded or not part of the training set
					addValidationRecord(record,r);
			}else
				if(exclude && inValidation.contains(rid)) // we have to predict it as part of this validation set...
					addValidationRecord(record,r);

		}


		client.out.println("eventually: training  " + trainingSetSize + " valid " + validationRowNums.size());
	}

	private void addTrainingRecord(int record, int repetitions, WorkflowNodeData wndOriginal)
	{
		// all duplicates added together
		for (int l = 0; l < 2; l++)
			// we have two ports: one with descriptors and one with values; conditions became part of data on the selection step
			for (int i = 0; i < repetitions; i++)
				wndInput.ports.get(l).addRow(wndOriginal.ports.get(l).getRow(record));
		if (repetitions > 1)
			trainingRowRepetitions.put(record, repetitions);
	}

	private void addValidationRecord(int record, AbstractDataRow r)
	{
		validationRowNums.add(record); // this is record that will be predicted
		wndInput.ports.get(0).addRow(r);
		if(wndInput.ports.size() == 1) wndInput.ports.add(new DataTable()); // FIX to have a possibility of a normal apply of models
		wndInput.ports.get(1).addStubRow(); // no values, of course
	}

}
