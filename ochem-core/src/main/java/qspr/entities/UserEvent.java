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

import java.sql.Timestamp;
import java.text.SimpleDateFormat;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

@Entity
@XmlRootElement(name = "user-event")
public class UserEvent
{
	@Id
	@Column(name = "ue_id")
	@GeneratedValue(strategy = GenerationType.AUTO)
	@XmlAttribute
	public Long id;
	
	@ManyToOne
	@JoinColumn(name = "session_id")
	@XmlTransient
	public Session session;

	@Column
	public String type;
	
	@Column
	public String comment;
	
	@Column(name = "xml_description")
	public String descriptionXml;
	
	@Column
	@XmlTransient
	public Timestamp time;
	
	@XmlElement
	protected String getUsername()
	{
		if (session.user == null)
			return "Guest (" + session.ipAddress + ")";
		else
			return session.user.login;
	}
	
	@XmlElement
	protected String getDate()
	{
		return new SimpleDateFormat("EEEE, dd.MM.yyyy").format(time);
	}
	
	@XmlElement
	protected String getTime()
	{
		return new SimpleDateFormat("HH:mm").format(time);
	}
	
	
}
