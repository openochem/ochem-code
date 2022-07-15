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

package com.eadmet.business;

import java.io.Serializable;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Molecule;
import qspr.entities.PendingTask.TaskType;
import qspr.export.ExportableMolecule;
import qspr.export.ExportableSetConfiguration;
import qspr.metaserver.configurations.CorinaConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.StandartizationOptions;
import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.metaserver.protocol.Task;
import qspr.modelling.AbstractTaskProcessor;
import qspr.modelling.CDSModelProcessor;
import qspr.modelling.ModelApplierAttachment;
import qspr.modelling.ModelProcessor;
import qspr.util.ExportThread;
import qspr.util.Operation;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.NodeConfiguration;
import qspr.workflow.datatypes.WorkflowConfiguration;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

public class DescriptorsCalculatorProcessor extends AbstractTaskProcessor
{
	public StandartizationOptions structureStandartisation = new StandartizationOptions();
	public StructureOptimisationConfiguration structureOptimisation;
	public DescriptorsConfiguration descConfig;
	public Basket basket;

	public DataTable dtDescriptors;

	private DescriptorsCalculatorAttachment attachment;

	public DescriptorsCalculatorProcessor() {
		setOperationID(Operation.generateID());
		taskClass = TaskType.DESCRIPTOR_CALCULATION;
	}

	@Override
	protected String getTaskType()
	{
		return QSPRConstants.Workflow;
	}

	@Override
	protected Serializable getTaskConfiguration() throws Exception
	{
		/**
		 * This enforces corina optimization for 3D descriptors.
		 * If corina is not available, we have a problem. 
		 * No optimization selected -> corina optimization -> not available -> no model
		 * on the other hand, e.g. CDK needs 3D structures else throws exception
		 */ 
		if (descConfig.requires3D() && structureOptimisation == null)
			structureOptimisation = new CorinaConfiguration();

		NodeConfiguration selectConfigStud = new NodeConfiguration("descriptors-selector", null);
		selectConfigStud.skipNode = true;

		WorkflowConfiguration wc = new WorkflowConfiguration("CDS")
		.addNodeConfiguration(new NodeConfiguration("descriptors-processor", descConfig))
		.addNodeConfiguration(selectConfigStud);
		CDSModelProcessor.addMoleculeProcessingConfiguration(wc.nodesConfiguration, structureStandartisation, structureOptimisation);

		String descTypes = descConfig.types.toString();
		taskName = descTypes.length() > 255 ? descTypes.substring(0, 255) : descTypes;

		return wc;
	}

	@Override
	protected Serializable getTaskData() throws Exception
	{
		ModelApplierAttachment.fetchMoleculesFromDatabase(basket);
		DataTable mol = ModelProcessor.basketToSDFTable(basket);
		if(mol.getRowsNoErrorsSize() == 0) throw new UserFriendlyException("No molecules were provided");
		return new WorkflowNodeData(mol);
	}

	@Override
	public DescriptorsCalculatorAttachment getAttachment()
	{
		if (attachment != null)
			return attachment;

		attachment = new DescriptorsCalculatorAttachment();
		attachment.setWorkData(basket);
		attachment.descConfig = descConfig;
		attachment.structureOptimisation = structureOptimisation;
		attachment.structureStandartization = structureStandartisation;

		return attachment;
	}

	@Override
	protected void onTaskReceived(Task task) throws Exception
	{
		dtDescriptors = WorkflowNodeData.fromTask(task).ports.get(0);
	}

	@Override
	protected void restoreFromAttachment(Serializable attachment)
	{
		DescriptorsCalculatorAttachment a = (DescriptorsCalculatorAttachment) attachment;
		descConfig = a.descConfig;
		structureOptimisation = a.structureOptimisation;
		structureStandartisation = a.structureStandartization;
		basket = a.getWorkData(false);
	}

	public ExportThread getExportThread(String format, ExportableSetConfiguration expConfig)
	{
		final DataTable dtDescriptors = this.dtDescriptors;
		return new ExportThread(format, expConfig)
		{
			@Override
			public void generateData() throws Exception
			{
				eData.setDescriptors(dtDescriptors);

				dtDescriptors.reset();
				while (dtDescriptors.nextRow())
				{
					if (eData.exportableMolecules.size() % 1000 == 0 && eData.exportableMolecules.size() > 0)
						Globals.restartAllTransactions(true);
					if (eData.exportableMolecules.size() % 100 == 0)
						setStatus("Preparing item " + eData.exportableMolecules.size());
					ExportableMolecule eMol = new ExportableMolecule();
					eData.addMolecule(eMol);

					ExperimentalProperty ep = basket.entries.get(dtDescriptors.currentRow).ep;

					if (ep.id != null)
					{
						ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, ep.id);
						eMol.setExperimentalProperty(ep);
					}
					else
					{
						Molecule mol = Repository.molecule.getMolecule(ep.molecule.id);
						eMol.setMolecule(mol);
					}

					eMol.descriptors = dtDescriptors.getCurrentRow();
				}

				//eData.supplementaryData.put("Configuration", model.configurationXml);
				setFileName("descriptors");
			}
		};
	}

	@Override
	public int bonusMultiplier() {
		return QSPRConstants.DESCRIPTORS_BONUS;
	}
}

class DescriptorsCalculatorAttachment extends ModelApplierAttachment
{
	StandartizationOptions structureStandartization = new StandartizationOptions();
	StructureOptimisationConfiguration structureOptimisation;
	DescriptorsConfiguration descConfig;

	private static final long serialVersionUID = 1L;
}
