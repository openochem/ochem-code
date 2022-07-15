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

import java.text.SimpleDateFormat;
import java.util.Calendar;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.exception.ConstraintViolationException;

import qspr.exception.MarshallableException;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;

@XmlRootElement(name = "message")
public class Alert 
{
	@XmlElement
	public String message = "";
	
	@XmlAttribute
	public String context;
	
	@XmlAttribute
	public String type;
	
	@XmlAttribute
	public String title = "";
	
	@XmlElement
	public String time;
	
	
	@XmlTransient
	public Object attachment;
	
	public Alert(String msg)
	{
		message = msg;
	}
	
	public Alert(String msg, String context)
	{
		message = msg;
		this.context = context;
	}
	
	public static Alert custom(String type, String title, String msg)
	{
		Alert alert = new Alert();
		alert.message = msg;
		alert.title = title;
		alert.type = type;
		return alert;		
	}
	
	public static Alert Debug(String _context, String msg)
	{
		Alert alert = new Alert();
		alert.context = _context;
		alert.message = msg;
		alert.title = msg;
		alert.type = "debug";
		
		return alert;
	}
	
	public static Alert Error(String _context, String msg)
	{
		Alert alert = Debug(_context, msg);
		alert.type = "error";
		
		return alert;
	}
	
	public static Alert Exception(Exception e)
	{
		Alert alert = new Alert();
		alert.title = e.getMessage();
		alert.attachment = e;
		
		if (e instanceof ConstraintViolationException)
		{
			qspr.exception.ConstraintViolationException cve = new qspr.exception.ConstraintViolationException((ConstraintViolationException)e);
			e = cve;
			alert.title = "Cannot complete your operation,\nsince the record you are trying to delete or modify is referenced by other records!\n\nThe detailed cause:\n" + cve.getUserMessage();
		}
		if (e instanceof UserFriendlyException)
			alert.context = "friendly";
		if (e instanceof MarshallableException)
			alert.attachment = ((MarshallableException)e).getObject();
		
		
		alert.message = OCHEMUtils.exceptionToString(e);
		alert.type = "exception";
		
		if (alert.title == null)
			alert.title = alert.message;
		
		return alert;
	}
	
	@XmlElement
	public Object getAttachment()
	{
		// Allow possibility to attach XML-marshallable objects to Alerts
		if (!(attachment instanceof Exception))
			return attachment;
		return null;
	}
	
	public Alert()
	{
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat("dd-MM-yyyy HH:mm:ss");
		time = sdf.format(cal.getTime());
	}
}
