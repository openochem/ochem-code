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

package qspr.util;

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Map;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.Globals;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;

// Universal container to track long-running server-side operations
// Complemented with client side Javascript long-operation.js

// Midnighter on Dec 20, 2011

@XmlRootElement
public class Operation 
{
	public static Map<Long, Operation> operations = new HashMap<Long, Operation>();
	
	public long operationId;
	
	@XmlTransient
	public StatusTracker statusTracker = new StatusTracker();
	
	public String successURL; // optional
	
	public String userLogin;
	
	private transient long timeStarted;
	
	public static long generateID()
	{
		return Math.round(Math.random() * 100000000);
	}
	
	public void cancel() //Cancel implementation for wrapper threads
	{
		WrapperThread t = WrapperThread.getByOperationID(operationId);
		if (t != null)
			t.cancelRequested = true;
	}
	
	@XmlElement
	public String getStatus()
	{
		return this.statusTracker.get();
	}
	
	public void setStatus(String status)
	{
		this.statusTracker.set(status);
	}
	
	public Operation(long operationId)
	{
		this.operationId = operationId;
		this.timeStarted = Calendar.getInstance().getTimeInMillis();
		Operation.registerOperation(this);
		//operations.put(operationId, this);
	}
	
	public static void registerOperation(Operation operation)
	{
		if (Globals.userSession() != null)
			operation.userLogin = Globals.userSession().user == null ? QSPRConstants.ANONYMOUS : Globals.userSession().user.login;
		if (operations.get(operation.operationId) != null)
			throw new UserFriendlyException("Trying to register operation " + operation.operationId + " second time. Please, check your code for consistency.");
		operations.put(operation.operationId, operation);
	}
	
	public Operation()
	{
		
	}
	
	public static Operation getOperation(long id)
	{
		return operations.get(id);
	}
	
	@XmlElement(name = "time-started")
	protected String getTimeStarted()
	{
		return new SimpleDateFormat("d MMMM, HH:mm:ss").format(timeStarted);
	}
}
