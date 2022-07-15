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

package com.eadmet.moderation;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;

import qspr.entities.Property;
import qspr.entities.User;

/**
 * A report about the data awaiting approval for a particular property
 * @author midnighter
 *
 */
@XmlRootElement(name = "awaiting-approval-data")
public class AwaitingApprovalData 
{
	/**
	 * The property related to this report
	 */
	public Property property;
	
	/**
	 * The list of users who published the data awaiting approval
	 */
	@XmlElementWrapper(name = "users")
	@XmlElement(name = "user")
	public List<User> users = new ArrayList<User>();
	
	public long approvedRecords;
	public long awaitingApprovalRecords;
	public long approvedModels;
	public long awaitingApprovalModels;
	
	public boolean equals(AwaitingApprovalData data)
	{
		return data.property == property;
	}
	
	public AwaitingApprovalData()
	{
		
	}
	
	public AwaitingApprovalData(Property property)
	{
		this.property = property;
	}
	
	public void addUser(User u)
	{
		if (!users.contains(u))
			users.add(u);
	}
	
	
}
