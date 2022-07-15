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

package qspr.business;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import qspr.entities.Basket;
import qspr.entities.PendingTask.TaskType;
import qspr.metaserver.configurations.CompareScaffoldsConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.StandartizationOptions;
import qspr.metaserver.protocol.Task;
import qspr.modelling.AbstractTaskProcessor;
import qspr.modelling.ModelApplierAttachment;
import qspr.modelling.ModelProcessor;
import qspr.util.Operation;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.NodeConfiguration;
import qspr.workflow.datatypes.NodesConfiguration;
import qspr.workflow.datatypes.WorkflowConfiguration;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

/**
 * A processor that manages the SetCompare tasks
 * @author midnighter
 *
 */
public class SetComparisonProcessor extends AbstractTaskProcessor
{
	//	private static transient final Logger logger = Logger.getLogger(SetComparisonProcessor.class);

	public Basket basket1;
	public Basket basket2;
	public DescriptorsConfiguration descriptors;
	public WorkflowNodeData wndResult;
	public StandartizationOptions standartization = new StandartizationOptions(true);

	public SetComparisonProcessor() {
		setOperationID(Operation.generateID());
		taskClass = TaskType.SET_COMPARISON;
	}

	public void onTaskReceived(Task task) throws IOException, ClassNotFoundException
	{
		wndResult = WorkflowNodeData.fromTask(task);
		//basket1 = basket2 = null; // Release memory
	}

	@Override
	protected Serializable getTaskConfiguration() throws Exception
	{
		CompareScaffoldsConfiguration compareConfiguration = new CompareScaffoldsConfiguration();
		NodesConfiguration nodesConfiguration = new NodesConfiguration();
		nodesConfiguration.nodes.add(new NodeConfiguration("mol-standardizer1", standartization));
		nodesConfiguration.nodes.add(new NodeConfiguration("mol-standardizer2", standartization));
		nodesConfiguration.nodes.add(new NodeConfiguration("descriptors1", descriptors));
		nodesConfiguration.nodes.add(new NodeConfiguration("descriptors2", descriptors));
		nodesConfiguration.nodes.add(new NodeConfiguration("compare-scaffolds", compareConfiguration));
		return new WorkflowConfiguration("set-comparison",  nodesConfiguration);
	}

	@Override
	protected Serializable getTaskData() throws Exception
	{
		WorkflowNodeData wnd = new WorkflowNodeData(ModelProcessor.basketToSDFTable(basket1));
		wnd.addPort(ModelProcessor.basketToSDFTable(basket2));

		return wnd;
	}

	private SetComparisonAttachment attachment;

	@Override
	public SetComparisonAttachment getAttachment()
	{
		if (attachment != null)
			return attachment;

		attachment = new SetComparisonAttachment();
		attachment.set1Attachment.setWorkData(basket1);
		attachment.set2Attachment.setWorkData(basket2);
		attachment.descriptors = descriptors;

		return attachment;
	}

	@Override
	public void restoreFromAttachment(Serializable attachment)
	{
		// Do we need the molecules data now?
		boolean needMolecules = isRunning();

		SetComparisonAttachment scAttachment = (SetComparisonAttachment) attachment;
		basket1 = scAttachment.set1Attachment.getWorkData(needMolecules);
		basket2 = scAttachment.set2Attachment.getWorkData(needMolecules);
		descriptors = scAttachment.descriptors;
	}

	@Override
	protected String getTaskType()
	{
		return QSPRConstants.Workflow;
	}

	// Cache of descriptor-filtered molecule indices
	private Map<String, List<Integer>>[] molIndices = new Map[]{new HashMap(), new HashMap()};

	/**
	 * Return the indices of the molecules having a particular feature
	 * @param descriptor - the descriptor (feature) name
	 * @param set - the number of the set: 0 or 1
	 * @return
	 */
	public List<Integer> getMoleculeIndices(String descriptor, int set)
	{
		if (molIndices[set].get(descriptor) == null)
		{
			List<Integer> indices = new ArrayList<Integer>();
			DataTable dtDescriptors = wndResult.ports.get(set + 1);
			if (dtDescriptors.containsColumn(descriptor))
			{
				int colIndex = dtDescriptors.indexOfColumn(descriptor);
				dtDescriptors.reset();
				while (dtDescriptors.nextRow())
					if (dtDescriptors.getValue(colIndex) != null && Double.valueOf("" + dtDescriptors.getValue(colIndex)) != 0.0)
						indices.add(dtDescriptors.currentRow);

			}
			molIndices[set].put(descriptor, indices);
		}

		return molIndices[set].get(descriptor);
	}

	/**
	 * Get indices of the failed molecules
	 * @param set
	 * @return
	 */
	public List<Integer> getErrorIndices(int set)
	{
		List<Integer> indices = new ArrayList<Integer>();
		DataTable dtDescriptors = wndResult.ports.get(set + 1);
		dtDescriptors.reset();
		while (dtDescriptors.nextRow())
			if (dtDescriptors.getCurrentRow().isError())
				indices.add(dtDescriptors.currentRow);

		return indices;
	}

	/**
	 * Get the error message for a particular molecule
	 * @param set
	 * @param molNum
	 * @return
	 */
	public String getErrorMessage(int set, int molNum)
	{
		return wndResult.ports.get(set + 1).getRow(molNum).detailedStatus;
	}
}

class SetComparisonAttachment implements Serializable
{
	private static final long serialVersionUID = 1L;

	/**
	 * Molecules in the first set
	 */
	ModelApplierAttachment set1Attachment = new ModelApplierAttachment();

	/**
	 * Molecules in the second set
	 */
	ModelApplierAttachment set2Attachment = new ModelApplierAttachment();

	DescriptorsConfiguration descriptors;
}



