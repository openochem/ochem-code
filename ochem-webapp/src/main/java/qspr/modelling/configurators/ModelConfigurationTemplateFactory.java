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

package qspr.modelling.configurators;

import java.io.StringWriter;
import java.util.List;

import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;

import com.eadmet.utils.MAILERConstants;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Attachment;
import qspr.entities.Attachment.AttachmentType;
import qspr.entities.AttachmentSource;
import qspr.entities.Model;
import qspr.entities.ModelAttachment;
import qspr.entities.ModelConfigurationTemplate;
import qspr.entities.ModelConfigurationTemplate.TemplateType;
import qspr.entities.Session;
import qspr.export.ExportableModel;
import qspr.metaserver.configurations.BaggingConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.SelectionConfiguration;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.util.WrapperThread;

public class ModelConfigurationTemplateFactory 
{
	public static ModelConfigurationTemplate create(ExportableModel eModel, TemplateType type)
	{
		ModelConfigurationTemplate mcTemplate = new ModelConfigurationTemplate();
		eModel.trainingSetId = null;
		eModel.properties.clear();

		mcTemplate.modelTemplate = Repository.modelTemplate.getByName(eModel.method);

		mcTemplate.type = type;

		if (type == TemplateType.MODEL)
		{
			mcTemplate.descriptors = ((CDSConfiguration)eModel.attachment.configuration).descriptors.types.toString();
			mcTemplate.validation = eModel.attachment.protocol.toString();
		}

		ExportableModel eModelCopy = getPartialTemplate(eModel, type);
		mcTemplate.configuration = new Attachment<ExportableModel>(eModelCopy, AttachmentType.MARSHALABLE, AttachmentSource.ModelConfTemplate);
		mcTemplate.md5 = mcTemplate.configuration.getAttachmentReference();

		try
		{
			mcTemplate.configurationXml = getXml(eModelCopy);
		}
		catch (JAXBException e)
		{
			e.printStackTrace();
			mcTemplate.configurationXml = "Invalid XML";
		}
		mcTemplate.name = getName(mcTemplate);

		mcTemplate.session = Globals.userSession();
		return mcTemplate;
	}

	public static ModelConfigurationTemplate create(Model m)
	{
		return create(ExportableModel.create(m), TemplateType.MODEL);
	}

	public static ModelConfigurationTemplate create(Model m, TemplateType type)
	{
		return create(ExportableModel.create(m), type);
	}

	public static String getName(ModelConfigurationTemplate mcTemplate)
	{
		ModelAttachment attachment = ((ExportableModel)mcTemplate.configuration.getObject()).attachment;
		switch (mcTemplate.type)
		{
		case MODEL:
			String validation = "No Validation";
			if (attachment.protocol.validationConfiguration !=null && attachment.protocol.validationConfiguration instanceof BaggingConfiguration)
				validation = "Bagging" + attachment.protocol.validationConfiguration.ensembleSize;
			else if (attachment.protocol.validationConfiguration !=null)
				validation = "" + attachment.protocol.validationConfiguration.ensembleSize + "FCV";
			return mcTemplate.modelTemplate.name + "_" + mcTemplate.descriptors + "_" + validation;
		case METHOD:
			return mcTemplate.modelTemplate.name;
		case DESCRIPTORS:
			return ((CDSConfiguration)attachment.configuration).descriptors.types.toString();
		case VALIDATION:
			validation = "No Validation";
			if (attachment.protocol.validationConfiguration !=null && attachment.protocol.validationConfiguration instanceof BaggingConfiguration)
				validation = "Bagging" + attachment.protocol.validationConfiguration.ensembleSize;
			else if (attachment.protocol.validationConfiguration !=null)
				validation = "" + attachment.protocol.validationConfiguration.ensembleSize + "FCV";
			return validation;
		case SELECTION:
			return ((CDSConfiguration)attachment.configuration).selection.toString();
		default:
			throw new RuntimeException("Unknown template type");
		}
	}

	public static ExportableModel getPartialTemplate(ExportableModel eModel, TemplateType type) 
	{
		ExportableModel eModelCopy;
		if (type == TemplateType.MODEL)
		{
			eModelCopy = eModel;
			//mcTemplate.descriptors = ((CDSConfiguration)eModel.attachment.configuration).descriptors.types.toString();
			//mcTemplate.validation = eModel.attachment.protocol.toString();
		} else if (type == TemplateType.METHOD)
		{
			eModelCopy = new ExportableModel();
			eModelCopy.method = eModel.method;
			CDSConfiguration cdsConf = new CDSConfiguration();

			//if (!(eModel.attachment.configuration instanceof CDSConfiguration))
			//	throw new UserFriendlyException("Cannot use this type of models as a template: " + eModel.attachment.configuration.getClass().getSimpleName());

			cdsConf.modelConfiguration = ((CDSConfiguration)eModel.attachment.configuration).modelConfiguration;
			cdsConf.descriptors = new DescriptorsConfiguration();;
			cdsConf.selection = null;
			eModelCopy.attachment = new ModelAttachment();
			eModelCopy.attachment.datahandling = null;
			eModelCopy.attachment.protocol = null;
			eModelCopy.attachment.standartization = null;

			eModelCopy.attachment.configuration = cdsConf;
		}
		else if (type == TemplateType.DESCRIPTORS)
		{
			eModelCopy = new ExportableModel();
			CDSConfiguration cdsConf = new CDSConfiguration();
			cdsConf.descriptors = ((CDSConfiguration)eModel.attachment.configuration).descriptors;
			List<String> descriptors = ((CDSConfiguration)eModel.attachment.configuration).selection.descriptorAsStrings();
			if (descriptors != null && !descriptors.isEmpty())
			{
				cdsConf.selection = new SelectionConfiguration();
				cdsConf.selection.storeDescriptorAsStrings(descriptors);
			}
			else
				cdsConf.selection = null;

			eModelCopy.attachment = new ModelAttachment();
			eModelCopy.attachment.configuration = cdsConf;
			eModelCopy.attachment.datahandling = null;
			eModelCopy.attachment.protocol = null;
			eModelCopy.attachment.standartization = null;
		}
		else if (type == TemplateType.VALIDATION)
		{
			eModelCopy = new ExportableModel();
			eModelCopy.attachment = new ModelAttachment();
			eModelCopy.attachment.protocol = eModel.attachment.protocol;
			eModelCopy.attachment.datahandling = null;
			eModelCopy.attachment.standartization = null;
		}
		else if (type == TemplateType.SELECTION)
		{
			eModelCopy = new ExportableModel();
			CDSConfiguration cdsConf = new CDSConfiguration();
			cdsConf.selection = ((CDSConfiguration)eModel.attachment.configuration).selection;
			cdsConf.selection.fixedDescriptors = null;
			cdsConf.descriptors = new DescriptorsConfiguration();
			eModelCopy.attachment = new ModelAttachment();
			eModelCopy.attachment.configuration = cdsConf;
			eModelCopy.attachment.datahandling = null;
			eModelCopy.attachment.protocol = null;
			eModelCopy.attachment.standartization = null;
		}
		else
			throw new RuntimeException("Unsupported template type " + type);

		return eModelCopy;
	}

	public static String getXml(ExportableModel eModel) throws JAXBException
	{
		StringWriter writer = new StringWriter();
		Marshaller marshaller = Globals.jaxbContext.createMarshaller();
		marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT,
				new Boolean(true));
		marshaller.marshal(eModel, writer);
		return writer.toString()
				// a patch: in case of feature nets, do not store the sub-models in the XML, otherwise it can become huge
				.replaceAll("nodesConfiguration[^\\uFFFF]*nodesConfiguration", "nodesConfiguration>(skipped)</nodesConfiguration") 
				.replace(" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"", "");
	}

	public static void main(String[] args) throws Exception {
		new WrapperThread() {

			@Override
			public void wrapped() throws Exception {

				ModelConfigurationTemplate template = create(Repository.model.getById(13069L));
				template.session = Session.getFirstSession(MAILERConstants.ADMIN);
				Globals.session().save(template);
			}
		}.run();
	}
}
