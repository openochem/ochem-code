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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.business.BatchUpoadBrowserFilter;
import qspr.business.BatchUpoadBrowserFilter.StatusFilter;

import com.eadmet.batchupload.main.RecordPreview.PreviewRecordStatus;
import com.eadmet.batchupload.main.RecordPreview.UploadRecordStatus;

@XmlRootElement
public class UploadPreview extends ArrayList<RecordPreview> 
{
	private static final long serialVersionUID = 1L;
	
	public String md5;
	public String remappingSchemaMD5;
	
	@XmlTransient
	public Map<String, RecordPreview> map = new HashMap<String, RecordPreview>();
	
	public UploadPreview()
	{
		
	}
	
	public List<RecordPreview> getRecords() 
	{
		return this;
	}
	
	private boolean matchesFilter(BatchUpoadBrowserFilter filter, RecordPreview rp)
	{
		if (filter == null)
			return true; 
		
		if (filter.lines.size() > 0)
			if (!filter.lines.contains(rp.getRow()))
				return false;
		
		
		if (StatusFilter.all.equals(filter.status))
			return true;
		
		if (StatusFilter.invalid.equals(filter.status))
			return (rp.status != PreviewRecordStatus.valid);	
		
		if (StatusFilter.valid.equals(filter.status))
			return (rp.status == PreviewRecordStatus.valid);

		if (StatusFilter.error.equals(filter.status))
			return (rp.status == PreviewRecordStatus.error);
		
		if (StatusFilter.fatal_error.equals(filter.status))
			return (rp.status == PreviewRecordStatus.fatal_error);
		
		if (StatusFilter.warning.equals(filter.status))
			return (rp.status == PreviewRecordStatus.warning);

		if (StatusFilter.duplicate_internal.equals(filter.status))
			return (rp.status == PreviewRecordStatus.duplicate_internal);

		if (StatusFilter.duplicate_external.equals(filter.status))
			return (rp.status == PreviewRecordStatus.duplicate_external);
		
		return true;
	}
	
	public List<RecordPreview> getRecords(BatchUpoadBrowserFilter filter)
	{
		List<RecordPreview> result = new ArrayList<RecordPreview>();
		for (RecordPreview rp : this)
			if (matchesFilter(filter, rp))
				result.add(rp);
		return result;
	}
	
	public PreviewSummary getSummary()
	{
		Map<PreviewRecordStatus, Map<UploadRecordStatus,PreviewSummaryItem>> m = new HashMap<PreviewRecordStatus, Map<UploadRecordStatus,PreviewSummaryItem>>();
		for (RecordPreview rp : this) 
		{
			if (m.get(rp.status) == null)
				m.put(rp.status, new HashMap<UploadRecordStatus,PreviewSummaryItem>());
			if (m.get(rp.status).get(rp.uploadStatus) == null)
				m.get(rp.status).put(rp.uploadStatus, new PreviewSummaryItem(rp.status, rp.uploadStatus, 0L));
			m.get(rp.status).get(rp.uploadStatus).count++;
		}
		PreviewSummary l = new PreviewSummary();
		for (PreviewRecordStatus ps : m.keySet())
			for (UploadRecordStatus us : m.get(ps).keySet())
				l.add(m.get(ps).get(us));
		l.total = this.size();
		return l;
	}
}


