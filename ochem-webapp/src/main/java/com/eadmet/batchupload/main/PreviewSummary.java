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

package com.eadmet.batchupload.main;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import com.eadmet.batchupload.main.RecordPreview.PreviewRecordStatus;
import com.eadmet.batchupload.main.RecordPreview.UploadRecordStatus;

@XmlRootElement
public class PreviewSummary extends ArrayList<PreviewSummaryItem> 
{
	private static final long serialVersionUID = 1L;
	
	@XmlAttribute
	Integer total;

	@XmlElement
	public List<PreviewSummaryItem> getSummary() 
	{
		return this;
	}
	
}

@XmlRootElement
class PreviewSummaryItem
{
	@XmlAttribute
	PreviewRecordStatus status;
	@XmlAttribute
	UploadRecordStatus uploadStatus;
	@XmlAttribute
	Long count;
	
	public PreviewSummaryItem()
	{
		
	}
	
	public PreviewSummaryItem(PreviewRecordStatus status, UploadRecordStatus uploadStatus, Long count)
	{
		this.status = status;
		this.uploadStatus = uploadStatus;
		this.count = count;
	}
	
	public String toString()
	{
		String r = "("+status;
		if (uploadStatus != null)
			r += ", "+uploadStatus;
		r +=") = "+count;
		return r;
	}
}