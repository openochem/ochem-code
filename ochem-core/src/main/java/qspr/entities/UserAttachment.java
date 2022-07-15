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

import java.io.Serializable;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Attachment.AttachmentType;

/**
 * An entity to store user-specific attachments associated with a particular key (e.g., some user settings, dialog settings, etc)
 * @author midnighter
 */
@Entity	
public class UserAttachment 
{
	@Id
	@Column(name = "ua_id")
	@GeneratedValue
	@XmlAttribute
	public Long id;
	
	@ManyToOne
	@JoinColumn(name="user_id")
	@XmlElement(name="user")
	public User user;
	
	@ManyToOne(cascade={CascadeType.ALL}, fetch = FetchType.LAZY)
	@JoinColumn(name = "attachment_id")
	@XmlTransient
	public Attachment<Serializable> attachment;
	
	@Column(name = "attachment_key")
	public String key;

	/**
	 * Get an object for a user by a key
	 */
	public static Serializable get(User user, String key)
	{
		if (user == null)
			return null;
		
		UserAttachment ua = getUA(user, key);
		if (ua == null)
			return null;
		
		try
		{
			return ua.attachment.getObject();
		}
		catch (Exception e)
		{
			// If by some reason we can deserialize/unmarshall the object, treat it as if its not there
			// This will resolve the compatibility issues
			logger.warn("Could not unserialize the user attachment " + key, e);
			return null;
		}
	}
	
	/**
	 * Set an object for a user by a key
	 */
	public static void set(User user, String key, Serializable object, AttachmentType attachmentType)
	{
		if (user == null)
			return;
		
		UserAttachment ua = getUA(user, key);
		if (ua == null)
		{
			ua = new UserAttachment();
			ua.user = user;
			ua.key = key;
			ua.attachment = new Attachment<Serializable>(object, attachmentType, AttachmentSource.UserAttachment);
			Globals.session().save(ua);
		}
		else
		{
			ua.attachment = new Attachment<Serializable>(object, attachmentType, AttachmentSource.UserAttachment);
			Globals.session().saveOrUpdate(ua.attachment);
			Globals.session().saveOrUpdate(ua);
		}
	}
	
	private static UserAttachment getUA(User user, String key)
	{
		return (UserAttachment) Globals.session().createCriteria(UserAttachment.class)
				.add(Restrictions.eq("user", user))
				.add(Restrictions.eq("key", key))
				.uniqueResult();
	}
	
	private static final Logger logger = LogManager.getLogger(UserAttachment.class);

}
