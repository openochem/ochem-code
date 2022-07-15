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

package com.eadmet.batchupload.entityschema;

import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import com.eadmet.batchupload.main.BatchUploadMessage;
import com.eadmet.batchupload.main.BatchUploadMessage.BatchUploadMessageType;
import com.eadmet.batchupload.main.RecordStub.ColumnValue;

@XmlRootElement
public class RemappedValue 
{
	@XmlAttribute
	public String originalName;
	
	@XmlAttribute
	public String name;
	
	@XmlElement
	public Double minValue;

	@XmlElement
	public Double maxValue;
	
	@XmlElement
	public List<BatchUploadMessage> messages = new ArrayList<BatchUploadMessage>();

	public void setNames(String name)
	{
		this.originalName = this.name = name;
	}
	
	public boolean valid(boolean warningIsInvalid)
	{
		for (BatchUploadMessage m : messages)
			if ((m.type == BatchUploadMessageType.error) || (warningIsInvalid && (m.type == BatchUploadMessageType.warning)))
				return false;
		return true;
	}
	
	public String toString()
	{
		return name;
	}

	public void updateWith(ColumnValue value) 
	{
		value.value = name;
	}
}
