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

import java.util.ArrayList;
import java.util.List;

import qspr.entities.DataHandlingOptions;
import qspr.entities.ModelAttachment;
import qspr.entities.ModelConfigurationTemplate;
import qspr.entities.ModelConfigurationTemplate.TemplateType;
import qspr.export.ExportableModel;
import qspr.metaserver.configurations.DescriptorsConfiguration.MixturesProcessing;
import qspr.metaserver.configurations.ModelAbstractConfiguration.ScalingType;
import qspr.metaserver.configurations.DescriptorType;
import qspr.metaserver.configurations.DescriptorsEmptyConfiguration;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.configurations.NoDescriptors;
import qspr.metaserver.configurations.SelectionConfiguration;
import qspr.metaserver.configurations.StandartizationOptions;
import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.metaserver.configurations.ValidationConfiguration;
import qspr.metaserver.configurations.ValidationConfiguration.MixtureValidation;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.modelling.configurations.ExternalCondition;
import qspr.util.ClassCompressor;

/**
 * Generator of model templates based from a dot product of descriptors, models, selection and validation templates
 * It generates the templates by "crossing over" the four template types
 * 
 * This is one of the core classes for Comprehensive Modeling utility
 * @author midnighter
 */

public class CrossModelGenerator 
{
	public Long trainingSetId; 
	public List<Long> validationSetIds;
	public List<Integer> templateIds = new ArrayList<Integer>();
	public List<Long> propertyIds = new ArrayList<Long>();
	public List<Long> unitIds = new ArrayList<Long>();

	public boolean skipDublicates;
	public StandartizationOptions standartization = new StandartizationOptions();
	public DataHandlingOptions dataHandling = new DataHandlingOptions();
	public MixturesProcessing mixturesProcessing;
	public Boolean saveModels;
	public ScalingType scalingTypeX;
	public ScalingType scalingTypeY;
	public String additionaDescriptors;
	public String implicitValues;
	public Boolean sanitize;
	public Boolean crs;
	public boolean singleLearning;
	public List<ExternalCondition> externalDescriptors = null;
	public String versionOCHEM;
	public String taxonomy;
	public StructureOptimisationConfiguration optimisationConfiguration = null;
	public boolean featureNet;
	public MixtureValidation mixtureValidation = null;

	/**
	 * Generate a set of model templates based on the configured class fields
	 * @return
	 */
	public List<ExportableModel> generateModelTemplates()
	{
		List<ExportableModel> methods = new ArrayList<ExportableModel>();
		List<ExportableModel> descriptors = new ArrayList<ExportableModel>();
		List<ExportableModel> selections = new ArrayList<ExportableModel>();
		List<ExportableModel> validations = new ArrayList<ExportableModel>();
		List<ExportableModel> models = new ArrayList<ExportableModel>();

		for (Integer id : templateIds) {
			ModelConfigurationTemplate template = ModelConfigurationTemplate.getById(id);
			ExportableModel eModel = (ExportableModel) template.configuration.getObject();
			if (TemplateType.METHOD == template.type)
				methods.add(eModel);
			else if (TemplateType.DESCRIPTORS == template.type)
				descriptors.add(eModel);
			else if (TemplateType.SELECTION == template.type)
				selections.add(eModel);
			else if (TemplateType.VALIDATION == template.type)
				validations.add(eModel);
		}

		for (ExportableModel validation : validations) {
			for (ExportableModel selection : selections) {
				for (ExportableModel method : methods) {
					for (ExportableModel descriptor : descriptors) {
						ExportableModel eModel = crossOverTemplates(method, descriptor, selection, validation);
						CDSConfiguration cdsConf = (CDSConfiguration)eModel.attachment.configuration;
						cdsConf.optimisationConfiguration = optimisationConfiguration;
						if(!cdsConf.descriptors.requires3D()) // no optimization if not required
							cdsConf.optimisationConfiguration = null;
						eModel.attachment.standartization = standartization;
						eModel.attachment.datahandling = dataHandling;

						if(!cdsConf.isCompatibleDescriptorsAndMethod(cdsConf.modelConfiguration))continue;

						if(cdsConf.modelConfiguration.isSupportDescriptors())
							cdsConf.descriptors.mixtures = mixturesProcessing;
						else
							cdsConf.descriptors.mixtures =  null; //no descriptors and no mixture processing

						cdsConf.modelConfiguration.saveModels = saveModels;

						if (scalingTypeX != null && ! (cdsConf.modelConfiguration instanceof NoDescriptors))
							cdsConf.modelConfiguration.scaleTypeX = scalingTypeX;

						if (scalingTypeY != null)
							cdsConf.modelConfiguration.scaleTypeY = scalingTypeY;

						if(externalDescriptors != null && cdsConf.modelConfiguration.isSupportConditions())
							cdsConf.conditions = externalDescriptors;

						if(mixtureValidation != null && eModel.attachment.protocol.validationConfiguration != null) { // only if deSalt is not active
							eModel.attachment.protocol.validationConfiguration.mixtureValidation = mixtureValidation;
							if(cdsConf.hasDescriptors() && cdsConf.descriptors.mixtures == null 
									&& cdsConf.modelConfiguration.isSupportDescriptors()) // just to indicate to algorithm that we work with mixtures
								cdsConf.descriptors.mixtures = MixturesProcessing.FRACTION;
						}

						cdsConf.modelConfiguration.setVersion(versionOCHEM, true);

						models.add(eModel);
					}
				}
			}
		}

		return models;
	}

	public ExportableModel crossOverTemplates(ExportableModel emMethod, ExportableModel emDescriptors, ExportableModel emSelection, ExportableModel emValidation)
	{
		ExportableModel model = new ExportableModel();
		model.attachment = new ModelAttachment();
		CDSConfiguration cdsConfiguration = new CDSConfiguration();
		model.attachment.configuration = cdsConfiguration;

		model.method = emMethod.method;
		cdsConfiguration.modelConfiguration = ((CDSConfiguration)emMethod.attachment.configuration).modelConfiguration;
		cdsConfiguration.descriptors = ((CDSConfiguration)emDescriptors.attachment.configuration).descriptors;

		if(additionaDescriptors != null){
			boolean found = false;
			for(DescriptorType type: cdsConfiguration.descriptors.types)
				if(additionaDescriptors.equals(type.type))found = true;
			if(!found){
				DescriptorType type = new DescriptorType(additionaDescriptors,new DescriptorsEmptyConfiguration());
				type.markUncachedAsErrors = true;
				cdsConfiguration.descriptors.types.add(type);
				cdsConfiguration.descriptors.allowMerge = true;
			}
		}

		cdsConfiguration.selection = ((CDSConfiguration)emSelection.attachment.configuration).selection;

		// Include the fixed descriptors, if provided in the descriptors template
		SelectionConfiguration descriptorsSelection = ((CDSConfiguration)emDescriptors.attachment.configuration).selection;
		if (descriptorsSelection.getDescriptorsSize() > 0)
		{
			cdsConfiguration.selection = (SelectionConfiguration) ClassCompressor.cloneObject(cdsConfiguration.selection); // clone the template to preserve it
			cdsConfiguration.selection.storeDescriptorAsStrings(descriptorsSelection.descriptorAsStrings());
		}

		if(cdsConfiguration.modelConfiguration.isSupportConditions()) // conditions are included
			cdsConfiguration.conditions = ((CDSConfiguration)emDescriptors.attachment.configuration).conditions;

		model.attachment.protocol = emValidation.attachment.protocol;

		model.attachment.datahandling = emMethod.attachment.datahandling;

		if(implicitValues != null)
			model.implicitValues = implicitValues;

		if(taxonomy != null)
			model.taxonomy = taxonomy;

		if(cdsConfiguration.modelConfiguration instanceof MultiLearningAbstractConfiguration && singleLearning) {
			((MultiLearningAbstractConfiguration)cdsConfiguration.modelConfiguration).noMultiLearning = true;
		}

		if(cdsConfiguration.modelConfiguration instanceof NoDescriptors) {
			if(sanitize)
				((NoDescriptors)cdsConfiguration.modelConfiguration).setSanitize();
		}

		if(featureNet) {
			cdsConfiguration.modelConfiguration.setFeatureNet(true);

			if(cdsConfiguration.modelConfiguration instanceof ValidationConfiguration) {
				ValidationConfiguration valid = (ValidationConfiguration) cdsConfiguration.modelConfiguration;
				valid.taskConfiguration.setFeatureNet(true);

				if(valid.taskConfiguration instanceof MultiLearningAbstractConfiguration)
					((MultiLearningAbstractConfiguration)valid.taskConfiguration).noMultiLearning = true;

			}

		}

		return model;
	}
}
