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

package com.eadmet.useractions;

import java.io.StringWriter;
import java.sql.Timestamp;
import java.util.Calendar;

import qspr.Globals;
import qspr.entities.UserEvent;

/**
 * A factory for the UserEvent class
 * @author midnighter
 *
 */
public class EventFactory
{
	/**
	 * Document a particular user action
	 * @param type type of the action
	 * @param actionData an optional object containing the details of the action
	 * @param comment an optional textual comment of the action. Can also be generated automatically from @param actionData
	 */
	public static void document(String type, AbstractUserAction actionData, String comment)
	{
		try
		{
			if (Globals.userSession() == null || (Globals.userSession().user != null && Globals.userSession().user.login.startsWith("Test")))
				return;
			
			UserEvent event = new UserEvent();
			event.type = type;
			event.session = Globals.userSession();
			event.time = new Timestamp(Calendar.getInstance().getTimeInMillis());
			event.comment = comment;
			
			if (actionData != null)
			{
				StringWriter writer = new StringWriter();
				Globals.jaxbContext.createMarshaller().marshal(actionData, writer);
				event.descriptionXml = writer.toString();
				event.comment = actionData.getLogLine();
			}
		
			Globals.session().save(event);
		}
		catch (Exception e)
		{
			// Never fail because of this. A failure in documentation should not be blocking.
			e.printStackTrace();
		}
	}
	
	public static void document(String type, AbstractUserAction actionData)
	{
		document(type, actionData, null);
	}
}
