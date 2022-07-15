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
import java.util.Calendar;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.annotations.Loggable;
import qspr.util.DynaWrap;

@Entity
@XmlRootElement(name = "message")
public class Message 
{
	public static final Boolean Unread_Message = false;
	public static final Boolean Read_Message = true;
	
	@Id
	@Column(name="message_id")
	@GeneratedValue
	@XmlAttribute
	public Long id;
	
	@Column
	public String subject;
	
	@Column
	public String body;
	
	@Column(name = "is_read")
	public boolean isRead;
	
	@ManyToOne
	@JoinColumn(name = "receiver_id")
	@XmlTransient
	@Loggable
	public User receiver;
	
	@ManyToOne
	@JoinColumn(name = "sender_id")
	@XmlTransient
	@Loggable
	public User sender;
	
	@ManyToOne
	@JoinColumn(name = "original_message_id")
	@XmlTransient
	public Message original_message;
	
	@Column(name="send_time")
	@XmlTransient
	public Timestamp time;
	
	@Column
	public String text;
	
	@XmlElement(name="sendername")
	public String getSenderName()
	{
		if (sender != null) {
			if (sender.isExtended()) {
				DynaWrap extended = sender.dynaWrapped();
				return extended.getString("firstName") + " " + extended.getString("lastName");
			} else {
				return sender.login;
			}
		}
		else
			return null;
	}
	
	@XmlElement(name="receivername")
	public String getReceiverName()
	{
		if (receiver != null) {
			if (receiver.isExtended()) {
				DynaWrap extended = receiver.dynaWrapped();
				return extended.getString("firstName") + " " + extended.getString("lastName");
			} else {
				return receiver.login;
			}
		}
		else
			return null;
	}
	
	@XmlElement(name = "sender")
	public String getSender()
	{
		if (sender != null)
			return sender.login;
		else
			return null;
	}
	
	@XmlElement(name = "receiver")
	public String getReceiver()
	{
		if (receiver != null)
			return receiver.login;
		else
			return null;
	}
	
	@XmlAttribute(name = "org_message_id")
	public Long getOriginalMessageid()
	{
		return original_message == null ? null : original_message.id;
	}
	
	@XmlElement
	public String getTime()
	{
		return (time != null)?ThreadScope.get().fullDateFormat.format(time):"";
	}
	
	@XmlAttribute(name = "login")
	public String getLoginName()
	{
		return Globals.userSession().user == null ? null : Globals.userSession().user.login;
	}
	
	public Message()
	{
		isRead = Unread_Message;
		time = new Timestamp(Calendar.getInstance().getTimeInMillis());
		original_message = this;
	}
}
