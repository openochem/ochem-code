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

import java.io.Serializable;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.entities.Attachment;
import qspr.entities.AttachmentSource;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.metaserver.configurations.BaggingConfiguration;
import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.metaserver.configurations.LabelWeighting;
import qspr.metaserver.configurations.LabelWeighting.PropertyWeighting;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.configurations.SelectionConfiguration;
import qspr.metaserver.configurations.StandartizationOptions;
import qspr.metaserver.configurations.ValidationConfiguration;
import qspr.workflow.utils.QSPRConstants;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.modelling.configurations.CDSModelData;
import qspr.modelling.configurations.ExternalCondition;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.NodeConfiguration;
import qspr.workflow.datatypes.NodesConfiguration;
import qspr.workflow.datatypes.WorkflowConfiguration;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.exceptions.UserFriendlyException;

public class CDSModelProcessor extends BasicModelProcessor
{
	private static transient final Logger logger = LogManager.getLogger(CDSModelProcessor.class);
	private static final int PROPERTY_OPTIONS = 200; // In reordering
	public boolean methodNeedsStructures = false;
	public int applierRepostSize = 0; // Allows to configure parallelization of predictions

	@Override
	public String getTaskType()
	{
		return QSPRConstants.Workflow;
	}

	@Override
	protected void prepare() throws Exception
	{
		super.prepare();
		if(getModelConfiguration().modelConfiguration instanceof MultiLearningAbstractConfiguration ) {

			MultiLearningAbstractConfiguration conf = (MultiLearningAbstractConfiguration)getModelConfiguration().modelConfiguration;

			if(model.implicitValues != null)
				conf.implicitValues = model.generateImplicitValues();
		}
	}

	private ModelAbstractConfiguration getModelData(WorkflowNodeData teacherResponse) throws Exception
	{
		if(teacherResponse.ports.get(1).getRowsSize()==0)return null;
		// Stub for black box models
		return (ModelAbstractConfiguration) teacherResponse.ports.get(1).getValue(0, 0);
	}

	private NodeConfiguration getApplierNodeConfiguration(Object modelData) throws Exception
	{
		// Stub for black box models

		if (modelData instanceof ModelAbstractConfiguration)
		{
			((ModelAbstractConfiguration)modelData).predictionScenario = predictionScenario;
			ModelAbstractConfiguration cfg = (ModelAbstractConfiguration)modelData;
			if (!cfg.isModelSaved() && !(cfg.versionOCHEM != null && cfg.versionOCHEM.equals("-999"))) // special to avoid Mock Model failure
				throw new UserFriendlyException("The model you selected cannot be applied to predict new compounds: the model was not saved");
		}

		if (model.attachment.getObject().protocol.validationConfiguration instanceof BaggingConfiguration)
		{
			NodeConfiguration nc = new NodeConfiguration("method-node", (Serializable) modelData);
			nc.taskType = "Bagging";
			return nc;
		}

		NodeConfiguration nc = new NodeConfiguration("method-node", (Serializable) modelData);
		nc.taskType = getTrainingMethodName(model.template.name);
		return nc;
	}

	private LabelWeighting reorderLabelWeighting(LabelWeighting lw)
	{
		if (lw == null)
			return null;

		if (lw.propertyWeights == null)
			return lw;

		LabelWeighting lwReordered = new LabelWeighting(lw);

		for (int i = 0; i < model.modelMappings.size(); i++)
		{
			Property p = model.modelMappings.get(i).property;
			PropertyWeighting pw = lw.getProperty(p.getName());

			if (pw == null)
				throw new UserFriendlyException("Missed label weighting in the model configuration: " + p.getName());

			PropertyWeighting pwNew = new PropertyWeighting(pw.name, pw.weight);
			lwReordered.addNewProperty(pwNew);

			PropertyOption[] options = new PropertyOption[PROPERTY_OPTIONS]; // ordered array of options
			int num = 0;
			for (Long optionId : model.attachment.getObject().optionsMapping.keySet())
			{
				Long mappedClassValue = model.getMappedOption(optionId);
				PropertyOption pOption = (PropertyOption) Globals.session().get(PropertyOption.class, optionId);
				if (pOption.property.equals(p))
				{
					options[mappedClassValue.intValue()] = pOption;
					num++;
				}
				if(num>=options.length)
					throw new UserFriendlyException("Increase number of property options " + num + " at reorderLabelWeighting");
			}

			for (int k = 0; k < num; k++)
				if (pw.getClass(options[k].name) == null)
					throw new UserFriendlyException("Incompatible weighting parameters in the model configuration: unknown class " + options[k].name);
				else
					lwReordered.addClass(p.getName(), options[k].name, pw.getClass(options[k].name).weight);

			// Reorder the cost matrix properly
			if (lw.useCostMatrix())
			{
				int classCount = pwNew.classesWeights.size();
				for (int j = 0; j < classCount; j++)
					for (int k = 0; k < classCount; k++)
					{
						// Addres classes by theiur names, not positions
						String c1 = pwNew.classesWeights.get(j).name;
						String c2 = pwNew.classesWeights.get(k).name;
						pwNew.setCostMatrix(c1, c2, pw.getCostMatrix(c1, c2));
					}

			}
		}

		logger.info(lwReordered);

		return lwReordered;
	}

	@Override
	protected final void onTeacherFinished(WorkflowNodeData teacherResponse) throws Exception
	{
		CDSConfiguration cdsConfiguration = (CDSConfiguration) model.attachment.getObject().configuration;
		CDSModelData cdsData = new CDSModelData();

		dtTrainingSetDescriptors = teacherResponse.getPort("descriptors");
		//cdsData.chosenDescriptors = new StringList(dtTrainingSetDescriptors.columns);
		cdsData.selectionConfiguration = teacherResponse.getPort("selection-configuration") !=null ? (SelectionConfiguration) teacherResponse.getPort("selection-configuration").getValue(0, 0) : null;
		cdsData.descriptors = cdsConfiguration.descriptors;
		cdsData.optimisationConfiguration = cdsConfiguration.optimisationConfiguration;
		cdsData.methodSpecificData = getModelData(teacherResponse);

		cdsData.methodSpecificData.setVersion(cdsConfiguration.modelConfiguration.versionOCHEM, true);

		if (model.getModelData(false) == null)
			model.setModelData(cdsData, false); //Original

		model.setModelData(cdsData, true); //Recalculated

		model.fullEquation = ((CDSModelData)model.getModelData(true)).getFullEquation();

		if(dtTrainingSetDescriptors != null)
			model.calcDescriptors = new Attachment<DataTable>(dtTrainingSetDescriptors.compress(), AttachmentSource.Descriptors);

		super.onTeacherFinished(teacherResponse);
	}

	@Override
	protected void onApplierFinished(WorkflowNodeData applierResponse) throws Exception
	{
		super.onApplierFinished(applierResponse);
		//dtValidationSetDescriptors = applierResponse.getPort("descriptors");
	}

	private NodesConfiguration getApplierNodesConfiguration(boolean recalculated, boolean forceRecalculateDescriptors) throws Exception 
	{
		CDSModelData cdsData = (CDSModelData) model.getModelData(recalculated);
		if (cdsData == null)
			throw new UserFriendlyException("The selected model does not contain model data and therefore cannot be applied to new compounds.\n");
		NodesConfiguration nodesConfiguration = new NodesConfiguration();

		// optimization
		if(forceRecalculateDescriptors){
			cdsData.descriptors.forceUpdateDescriptorCache = true;		
			if(cdsData.optimisationConfiguration != null)
				cdsData.optimisationConfiguration.bypassCache = true; // we also recalculate structures, i.e. all steps
		}

		addMoleculeProcessingConfiguration(nodesConfiguration, model.attachment.getObject().standartization, cdsData.optimisationConfiguration);
		nodesConfiguration.nodes.add(new NodeConfiguration("descriptors-processor", cdsData.descriptors));
		nodesConfiguration.nodes.add(new NodeConfiguration("descriptors-selector", cdsData.selectionConfiguration));

		if (cdsData.methodSpecificData == null)
			throw new UserFriendlyException(" This model was not saved and cannot be applied to new data");

		NodeConfiguration nodeConfiguration = getApplierNodeConfiguration(cdsData.methodSpecificData);
		if (nodeConfiguration != null)
			nodesConfiguration.nodes.add(nodeConfiguration);

		return nodesConfiguration;
	}

	public static void addMoleculeProcessingConfiguration(NodesConfiguration nc, StandartizationOptions standartization, StructureOptimisationConfiguration optimization)
	{
		NodeConfiguration ncStandartizarion = new NodeConfiguration("mol-standardizer", standartization);
		ncStandartizarion.taskType = QSPRConstants.MolStandartizer;
		nc.nodes.add(ncStandartizarion);

		NodeConfiguration ncOptimisation = new NodeConfiguration("mol-optimiser", optimization);
		nc.nodes.add(ncOptimisation);

		if (StructureOptimisationConfiguration.isDefaultConfiguration(optimization))
		{
			ncOptimisation.skipNode = true;
			// A piece of legacy logic: if we do not optimize structures (Hydrogens will be added there), add hydrogens explicitly
			standartization.addExplicitHydrogensWith = standartization.getDefault();
			standartization.dearomatizeWith = standartization.getDefault();
		}
		else
			ncOptimisation.taskType = optimization.getTaskType();
	}

	private NodesConfiguration getTeacherNodesConfiguration() throws Exception 
	{
		CDSConfiguration configuration = (CDSConfiguration) model.attachment.getObject().configuration;
		NodesConfiguration nodesConfiguration = new NodesConfiguration();

		addMoleculeProcessingConfiguration(nodesConfiguration, model.attachment.getObject().standartization, configuration.optimisationConfiguration);
		nodesConfiguration.nodes.add(new NodeConfiguration("descriptors-processor", configuration.descriptors));
		nodesConfiguration.nodes.add(new NodeConfiguration("descriptors-selector", configuration.selection));

		NodeConfiguration nodeConfiguration = getTeacherNodeConfiguration(configuration.modelConfiguration);
		if (nodeConfiguration != null)
			nodesConfiguration.nodes.add(nodeConfiguration);

		return nodesConfiguration;
	}

	private NodeConfiguration getTeacherNodeConfiguration(ModelAbstractConfiguration modelConfiguration)
	{
		if(modelConfiguration instanceof MultiLearningAbstractConfiguration)
			((MultiLearningAbstractConfiguration)modelConfiguration).labelWeighting = reorderLabelWeighting(((MultiLearningAbstractConfiguration)modelConfiguration).labelWeighting);
		modelConfiguration.setTheOptions(numOfOptions);	

		model.isConsistent();

		if (model.attachment.getObject().protocol.validationConfiguration != null)
		{
			// Bagging 
			ValidationConfiguration baggingConf = (ValidationConfiguration)model.attachment.getObject().protocol.validationConfiguration.getDeepCopy();

			baggingConf.taskConfiguration = modelConfiguration;
			baggingConf.taskName = getTrainingMethodName(model.template.name);

			NodeConfiguration nc = new NodeConfiguration("method-node", baggingConf);
			nc.taskType =  baggingConf instanceof BaggingConfiguration? "Bagging":"CrossValidation";

			return nc;
		}
		else
		{
			// No validation
			NodeConfiguration nc = new NodeConfiguration("method-node", modelConfiguration);
			nc.taskType = getTrainingMethodName(model.template.name);
			return nc;
		}	
	}

	@Override
	protected void preloadLazyData()
	{
		super.preloadLazyData();
		if (getModelConfiguration().hasConditions())
			for (ExternalCondition eDesc : getModelConfiguration().conditions)
			{
				eDesc.getProperty();
				eDesc.getUnit();
			}
	}

	private String getWorkflowName(boolean applier)
	{
		// kernel stuff are special case
		boolean kernel = model.template.name.equals("KPLS");
		String name = kernel ? model.template.name.toLowerCase() : "cds";
		name = name + (applier ? "-applier" : "-teacher");

		if (!kernel && methodNeedsStructures)
			name += "-with-structures";

		return name;
	}

	@Override
	CDSConfiguration getModelConfiguration()
	{
		return (CDSConfiguration) super.getModelConfiguration();
	}

	@Override
	public Task createApplierTask() throws Exception
	{
		return new Task(QSPRConstants.Workflow, getApplierConfiguration(true, false), null);
	}

	@Override
	protected Task createTeacherTask() throws Exception
	{
		return new Task(QSPRConstants.Workflow, getTeacherConfiguration(), null);
	}

	@Override
	public Serializable getApplierConfiguration(boolean recalculated, boolean forceRecalculateDescriptors) throws Exception
	{
		return new WorkflowConfiguration(getWorkflowName(true), getApplierNodesConfiguration(recalculated, forceRecalculateDescriptors)).setRepostSize(applierRepostSize);
	}

	@Override
	Serializable getTeacherConfiguration() throws Exception
	{
		return new WorkflowConfiguration(getWorkflowName(false), getTeacherNodesConfiguration());
	}


	private String getTrainingMethodName(String modelTemplateName)
	{
		if (modelTemplateName.equals("LibraryModel"))
			return "ASNN";
		else
			return model.template.name;
	}
}
