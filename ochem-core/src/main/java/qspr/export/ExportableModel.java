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

package qspr.export;

import java.io.ByteArrayOutputStream;
import java.io.Serializable;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import javax.xml.bind.JAXBException;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Hibernate;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.dao.Repository;
import qspr.entities.Attachment;
import qspr.entities.Attachment.AttachmentType;
import qspr.entities.AttachmentSource;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.ModelAttachment;
import qspr.entities.ModelMapping;
import qspr.entities.ModelMicroAttachment;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.util.ClassCompressor;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;

/**
 * Allows to marshall a model (configuration, training set name, property names) into an XMl file.
 * This class is not supposed to export the ready model itself but can be used to redevelop a model
 * 
 * @author midnighter
 *
 */
@XmlRootElement(name = "ochem-model")
public class ExportableModel 
{
	public String name;
	public String description;

	@XmlElementWrapper(name = "properties")
	@XmlElement(name = "property")
	public List<ModelProperty> properties = new ArrayList<ModelProperty>();

	public String method;

	public ModelAttachment attachment;

	public Long trainingSetId;
	public List<Long> validationSetId;


	public Long originalPublicId;
	public Long originalInternalId;
	public String exportDate;

	@XmlTransient
	public String implicitValues;
	public String taxonomy;
	
	public List<String> getPropertyOptions(int i) {
		return properties.get(i).options;
	}

	public static ExportableModel create(Model m)
	{
		if(m == null)throw new UserFriendlyException("Some of your models have been deleted. Refresh CM view before starting recalculation of models.");
		ExportableModel eModel = new ExportableModel();
		eModel.description = m.userDescription;
		eModel.originalPublicId = m.publicId;
		eModel.originalInternalId = m.id;
		eModel.name = m.featuredName != null ? m.featuredName : m.name;
		eModel.exportDate = new SimpleDateFormat("yyyy-MM-dd HH:mm").format(Calendar.getInstance().getTime());

		for (ModelMapping mm : m.modelMappings)
			eModel.properties.add(new ModelProperty(mm));
		eModel.method = m.template.name;
		eModel.trainingSetId = m.trainingSet.id;
		List<Basket> valSets = m.getValidationSets();
		if (valSets != null && !valSets.isEmpty())
		{
			eModel.validationSetId = new ArrayList<Long>();
			for (Basket valSet : valSets)
				eModel.validationSetId.add(valSet.id);
		}

		ModelAttachment attachment = m.attachment.getObject();
		eModel.attachment = new ModelAttachment();
		if (ClassCompressor.classLoader == null)
			ClassCompressor.classLoader = ThreadScope.get().threadClassLoader;
		eModel.attachment.configuration = ClassCompressor.cloneObject((Serializable)attachment.configuration);
		eModel.attachment.datahandling = attachment.datahandling;
		eModel.attachment.protocol = attachment.protocol;
		eModel.attachment.standartization = attachment.standartization;

		return eModel;
	}

	public Model createModel()
	{
		Model model = new Model();
		if (trainingSetId != null)
			model.trainingSet = (Basket) Globals.session().get(Basket.class, trainingSetId);
		if (validationSetId != null)
			for (Long id : validationSetId)
				model.addValidationSet((Basket) Globals.session().get(Basket.class, id));

		model.attachment = new Attachment<ModelAttachment>(attachment, AttachmentType.MARSHALABLE, AttachmentSource.Model);
		model.microattachment = new Attachment<ModelMicroAttachment>(new ModelMicroAttachment(), AttachmentType.MARSHALABLE, AttachmentSource.Model);

		model.template = Repository.modelTemplate.getByName(method);
		if (properties.isEmpty() && model.trainingSet != null)
			for (Property prop : model.trainingSet.getProperty())
				properties.add(new ModelProperty(prop, model.attachment.getObject()));

		for (ModelProperty mProperty : properties) 
		{
			ModelMapping mm = new ModelMapping();
			mm.property = Repository.property.getProperty(mProperty.name, false);
			if (mm.property == null)
				throw new UserFriendlyException("Cannot find property with name " + mProperty.name);
			Hibernate.initialize(mm.property.options);
			model.addModelMapping(mm);
			if (mm.property.isNumeric())
				if (mProperty.unit != null)
					mm.unit = Repository.unit.get(mProperty.unit, mm.property.unitCategory);
				else
					mm.unit = mm.property.defaultUnit;
			else
				mm.unit = mm.property.defaultUnit;
		}

		// this is a fix!!!!
		model.implicitValues = implicitValues;
		model.taxonomy = taxonomy;
		
		return model;
	}

	public String getMD5() throws JAXBException
	{
		ByteArrayOutputStream stream = new ByteArrayOutputStream();
		Globals.jaxbContext.createMarshaller().marshal(this, stream);
		return OCHEMUtils.getMD5(stream.toByteArray());
	}

	public String toString()
	{
		boolean hasMethod = attachment.configuration != null && ((CDSConfiguration)attachment.configuration).modelConfiguration != null;
		boolean hasDescriptors = attachment.configuration != null && ((CDSConfiguration)attachment.configuration).hasDescriptors();
		if (hasMethod && hasDescriptors)
			return attachment.toString();
		else if (hasMethod)
			return ((CDSConfiguration)attachment.configuration).modelConfiguration.toString();
		else if (hasDescriptors)
			return ((CDSConfiguration)attachment.configuration).descriptors.toString();
		else
			return attachment.toString();

	}
}

class ModelProperty
{
	@XmlAttribute
	public String name;
	@XmlAttribute
	public String unit;
	@XmlAttribute
	public Long id;

	@XmlElement(name = "option")
	public List<String> options;

	public ModelProperty()
	{

	}

	public ModelProperty(Property p, ModelAttachment a)
	{
		this.name = p.getName();
		this.unit = p.defaultUnit.getName();
		this.id = p.id;

		if (p.isQualitative())
		{
			options = new ArrayList<String>();
			for (Long key : a.optionsMapping.keySet())
			{
				PropertyOption po = PropertyOption.getById(key);
				if (po.property.id.equals(p.id))
				{
					while (options.size() <= a.optionsMapping.get(key))
						options.add(null);
					options.set(a.optionsMapping.get(key).intValue(), po.name);
				}
			}
		}
	}

	public void addOptions(ModelAttachment attachment) {

	}

	public ModelProperty(ModelMapping mm)
	{
		this(mm.property, mm.model.attachment.getObject());
	}
}
