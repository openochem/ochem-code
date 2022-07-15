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

import javax.xml.bind.annotation.XmlRootElement;

import com.eadmet.datacube.Metrics;

@XmlRootElement
public class DataMetrics implements Metrics
{
	public int approvedRecords;
	public int privateRecords;
	public int awaitingApproval;
	public int awaitingApprovalOverdue;
	public int rejected;
	public int awaitingModels;
	public int privateModels;
	public int approvedModels;
	
	@Override
	public void add(Metrics metrics)
	{
		DataMetrics operand = (DataMetrics) metrics;
		
		approvedRecords += operand.approvedRecords;
		privateRecords += operand.privateRecords;
		awaitingApproval += operand.awaitingApproval;
		awaitingApprovalOverdue += operand.awaitingApprovalOverdue;
		rejected += operand.rejected;
		awaitingModels += operand.awaitingModels;
		privateModels += operand.privateModels;
		approvedModels += operand.approvedModels;
	}
	
	public String toString() {
		return "" + awaitingApproval + " records awaiting approval, " + approvedModels + " approved models";
	}
}
