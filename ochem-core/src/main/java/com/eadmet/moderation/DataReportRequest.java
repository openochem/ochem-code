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

import javax.xml.bind.annotation.XmlRootElement;

import qspr.entities.Property;
import qspr.entities.User;

import com.eadmet.datacube.Subkonto;

/**
 * A request for a multileveled report on data records and models filtered by a particular criteria
 * 
 * @author midnighter
 *
 */
@XmlRootElement
public class DataReportRequest
{
	public boolean countUnapprovedRecords = true;
	public boolean countApprovedRecords;
	public boolean countRejectedRecords = true;
	public boolean countPrivateRecords;
	public boolean countModels = true;
	public Long basketId; //Limit the report only to a basket
	public Boolean original; //If not null, limits the report only to original/unoriginal data
	
	public User moderator;
	public User user;
	public List<Property> properties = new ArrayList<Property>();
	
	/**
	 * In case we want a grouping by a qualitative condition, this field specifies which exactly.
	 * Currently, only one condition-based grouping per query is supported.
	 * 
	 * Under construction.
	 */
	public Property qualitativeCondition;
	
	public List<Subkonto> groupingOrder;
}
