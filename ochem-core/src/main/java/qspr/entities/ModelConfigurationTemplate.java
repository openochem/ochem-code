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

package qspr.entities;

import java.io.StringWriter;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.EnumType;
import javax.persistence.Enumerated;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.Transient;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Criteria;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.util.AccessChecker;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;

@Entity
@XmlRootElement(name = "model-template")
public class ModelConfigurationTemplate 
{
	@Id
	@GeneratedValue
	@Column(name = "mct_id")
	@XmlAttribute
	public Integer id;
	
	@Column
	@XmlAttribute
	public String name;
	
	@ManyToOne
	@JoinColumn(name = "session_id")
	public Session session;
	
	@ManyToOne
	@JoinColumn(name = "model_template_id")
	public ModelTemplate modelTemplate;

	/**
	 * Used to avoid duplicated cfg for the same user
	 */
	public String md5;

	// ** UI purposes only ** //
	@XmlTransient
	public String descriptors;
	@XmlTransient
	public String validation;
	
	@XmlTransient
	@Transient
	public String configurationXml;

	@XmlTransient
	@ManyToOne(cascade={CascadeType.ALL}, fetch = FetchType.LAZY)
	@JoinColumn(name = "configuration")
	public Attachment configuration;
	
	@Column@XmlElement
	public boolean recommended;
	
	@ManyToOne
	@JoinColumn(name = "introducer_id")
	public User introducer;
	
	@Column(name = "public")
	public boolean isPublic;
	
	@Enumerated(EnumType.STRING)
	@Column
	public TemplateType type;
	
	public ModelConfigurationTemplate()
	{
		
	}
	
	public void requestModification()
	{
		if (isPublic)
			AccessChecker.requestSuperuserPrivileges();
		else
			if (introducer != null && !introducer.equals(Globals.userSession().user))
				throw new UserFriendlyException("You cannot modify this template");
	}
	
	public boolean hasConflicts()
	{
		Criteria criteria = Globals.session().createCriteria(ModelConfigurationTemplate.class).add(Restrictions.eq("md5", md5));
		if (id != null)
			criteria.add(Restrictions.ne("id", id));
		criteria.setProjection(Projections.count("id"));
		return ((Long)criteria.uniqueResult()) > 0;
	}
	
	public void checkConflicts()
	{
		if (hasConflicts())
			throw new UserFriendlyException("The template you are trying to save is a dublicate for another template");
	}
	
	public static ModelConfigurationTemplate getById(Integer id)
	{
		return (ModelConfigurationTemplate) Globals.session().get(ModelConfigurationTemplate.class, id);
	}
	
	@XmlElement(name = "full-xml")
	public String getFullXml() throws JAXBException
	{
		if (!Globals.getMarshallingOption(MarshallingOption.MODELTEMPLATE_FULLXML))
			return null;

		Marshaller marshaller = Globals.jaxbContext.createMarshaller();
		marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT,
				   new Boolean(true));
		StringWriter writer = new StringWriter();
		marshaller.marshal(configuration.getObject(), writer);
		
		return writer.toString();
	}
	
	public void updateHash()
	{
		if (isPublic || introducer == null)
			md5 = configuration.getAttachmentReference();
		else
			md5 = OCHEMUtils.getMD5(configuration.getAttachmentReference() + introducer.id);
	}
	
	@XmlElement
	public String getDescription()
	{
		try
		{
			return configuration.getObject().toString();
		}
		catch (Exception e)
		{
			return "Incompatible XML. API might have been updated and your template is no more valid";
		}
	}
	
	public static enum TemplateType {
		MODEL, 
		DESCRIPTORS, 
		METHOD,
		SELECTION,
		VALIDATION
	};
}
