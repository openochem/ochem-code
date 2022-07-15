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

package qspr.business;

import com.eadmet.exceptions.UserFriendlyException;

/**
 * A simple representation of the privileges set of a particular user to a particular object.
 * @author midnighter
 *
 */
public class Privileges
{
	public boolean canView;
	public boolean canEdit;
	
	public String entityType;
	
	public Privileges(String entityType)
	{
		this.entityType = entityType;
	}
	
	public void requestView()
	{
		if (!canView)
			throw new UserFriendlyException("You are not authorized to access this " + entityType);
	}
	
	public void requestModification()
	{
		if (!canEdit)
			throw new UserFriendlyException("You are not authorized to modify this " + entityType);
	}
}
