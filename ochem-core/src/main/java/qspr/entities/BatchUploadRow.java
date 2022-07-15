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

package qspr.entities;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.EnumType;
import javax.persistence.Enumerated;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

@Entity
@XmlRootElement(name = "batchuploadrow")
public class BatchUploadRow 
{
	@Id
	@Column(name = "batchuploadrow_id")
	@GeneratedValue
	@javax.persistence.OrderBy
	@XmlAttribute
	public Long id;
	
	@ManyToOne(fetch=FetchType.LAZY)
	@JoinColumn(name = "batchupload_id")
	@XmlTransient
	public BatchUpload batchUpload;
	
	@Column
	@XmlTransient
	public Integer sheet;
	
	@Column
	@XmlAttribute
	public String num;

	@Enumerated(EnumType.STRING)
	@Column
	@XmlAttribute
	public RowStatus status = RowStatus.undefined;
	

	@Column(name="detailed_status")
	@XmlElement
	public String detailedStatus;

	@Column
	@XmlTransient
	public String trace;
	
	@Enumerated(EnumType.STRING)
	@Column(name = "upload_action")
	@XmlAttribute
	public UploadAction action = UploadAction.undefined;
	
	@Enumerated(EnumType.STRING)
	@Column(name = "type")
	@XmlAttribute
	public RowType type = RowType.undefined;

	@ManyToOne
	@JoinColumn(name = "exp_property_id", nullable = true)
	@XmlTransient
	public ExperimentalProperty ep;

	@Column
	@XmlTransient
	public byte[] data;

	public enum RowStatus
	{
		undefined, preloaded, preprocessed, uploading, uploaded, skipped, merged, error;
		
		public String toReadableString()
		{
			switch (this)
			{
				case undefined: return "undefined";
				case preloaded: return "preloaded";
				case preprocessed: return "preprocessed";
				case uploading: return "uploading";
				case error: return "error";
				case skipped: return "skipped"; 
				case uploaded: return "uploaded";
				case merged: return "merged";
			}
			return "broken";
		}
		
		public String toString()
		{
			return toReadableString();
		}
	}
	
	public enum RowType
	{
		undefined, valid, internal_duplicate, external_duplicate, warning, error, updated;
		
		public String toReadableString()
		{
			switch (this)
			{
				case external_duplicate: return "duplicate";
				case internal_duplicate: return "internal duplicate";
				case warning: return "warning";
				case undefined: return "undefined";
				case error: return "error";
				case valid: return "valid";
				case updated: return "updated";
			}
			return "broken";
		}		
		
		public String toString()
		{
			return toReadableString();
		}
	}
	
	public enum UploadAction
	{
		save, skip, merge, update, undefined;
		public String toReadableString()
		{
			return toString();
		}
	}
	
	public void appendDetailedStatus(String status)
	{
		if (detailedStatus == null)
			detailedStatus = status;
		else
			detailedStatus += "\n"+status;
				
	}
}
