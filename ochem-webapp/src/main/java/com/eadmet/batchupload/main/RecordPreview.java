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
import javax.xml.bind.annotation.XmlTransient;

import qspr.entities.ExperimentalProperty;

import com.eadmet.batchupload.main.BatchUploadMessage.BatchUploadMessageType;

@XmlRootElement
public class RecordPreview implements Cloneable
{
	
	@XmlAttribute
	public long id = 0;
	
	
	public PreviewUploadAction action;
	public PreviewRecordStatus status;
	public UploadRecordStatus uploadStatus;
	
	
	@XmlElement
	public List<BatchUploadMessage> messages = new ArrayList<BatchUploadMessage>();
	
	@XmlElement
	public ExperimentalProperty ep; //Store it only temporarily, remove to conserve memory
	
	@XmlElement
	public ExperimentalProperty duplicate; //Store it only temporarily, remove to conserve memory
	
	
	//
	@XmlTransient
	int sheet = -1;
	
	@XmlTransient
	int row = -1;
	
	@XmlTransient
	int record = -1; // Coordinates of this "preview" in the "file"
	
	@XmlAttribute
	public int getRow() // People like 1-based numbering instead of 0-based
	{
		return row+1;
	}
	
	@XmlAttribute
	public int getRecord()
	{
		return record+1;
	}
	
	@XmlAttribute
	public int getSheet()
	{
		return sheet+1;
	}
	
	public enum PreviewUploadAction {skip, save, merge, put_original_to_basket};
	public enum PreviewRecordStatus {valid, fatal_error, error, warning, duplicate_external, duplicate_internal};
	public enum UploadRecordStatus {saved_valid, saved_error, saved_duplicate, fatal_error, merged, skipped, put_original_to_basket};
	
	
	public String getCoordinates()
	{
		return sheet+"_"+row+"_"+record;
	}
	
	public void filterMessages()
	{
		int i = 1;
		while (i < messages.size())
			if (messages.get(i).message.equals("Null pointer exception"))
				messages.remove(i);
			else
				i++;
	}
	
	public boolean equalByCoordinates(RecordPreview rp)
	{
		return getCoordinates().equals(rp.getCoordinates());
	}
	
	public void minimize()
	{
		ep = null;
		duplicate = null;
		if (messages.size() > 1)
			messages.subList(1, messages.size()).clear();
	}
	
	public int numErrorMessages()
	{
		int num = 0;
		for (BatchUploadMessage m : messages) 
			if (m.type == BatchUploadMessageType.error)
				num++;
		return num;
	}
	
	@Override
	public RecordPreview clone()
	{
		RecordPreview rp = new RecordPreview();
		rp.id = id;
		rp.action = action;
		rp.status = status;
		rp.sheet = sheet;
		rp.row = row;
		rp.record = record;
		return rp;
	}
}
