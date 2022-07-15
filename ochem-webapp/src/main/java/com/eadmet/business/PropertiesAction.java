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

package com.eadmet.business;

import java.util.ArrayList;
import java.util.List;

public class PropertiesAction
{
	public Long id;
	public Long parentId;
	public String name;
	public boolean isDirectory = false;
	public int propertyType = 0;
	public boolean isCondition = false;
	public Long unitCategoryId;
	public boolean confirmed = false;
	public String description;
	public String aliases;
	public Long defaultUnitId;
	public boolean isPublic = false;
	public boolean isApproved = false;
	public Double bonusPointsWeight;
	public List<Long> obligatoryConditionIds = new ArrayList<Long>();
	
	public PropertiesAction() 
	{
	}
}
