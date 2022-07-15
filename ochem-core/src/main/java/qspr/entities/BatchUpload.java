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

import java.sql.Timestamp;
import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.EnumType;
import javax.persistence.Enumerated;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.entities.BatchUploadRow.RowStatus;
import qspr.entities.BatchUploadRow.RowType;

@Entity
@XmlRootElement(name = "batchupload")
public class BatchUpload 
{
	@Id
	@Column(name = "batchupload_id")
	@GeneratedValue
	@javax.persistence.OrderBy
	@XmlAttribute
	public Long id;

	@ManyToOne
	@JoinColumn(name = "session_id")
	@XmlTransient
	public Session session;

	@Enumerated(EnumType.STRING)
	@Column
	@XmlAttribute
	public BatchUploadState state = BatchUploadState.init;

	@Column
	@XmlElement
	public String status;

	@Column(name="last_activity")
	@XmlTransient
	public Timestamp lastActivity;

	@ManyToOne
	@JoinColumn(name="file_attachment_id")
	@XmlTransient
	public Attachment file;

	@Column(name="file_name")
	@XmlElement
	public String fileName;

	@Column(name="uploadcontext")
	@XmlTransient
	public String context;

	@Column(name="uploadschema")
	@XmlTransient
	public String schema;

	@Column(name="uploadstructure_original")
	@XmlTransient
	public String originalStructure; //Move it to attachments...?

	@Column(name="uploadstructure_processed")
	@XmlTransient
	public String processedStructure;

	@XmlElement
	public String getLastActivity()
	{
		if (lastActivity == null)
			return "";
		Format formatter = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		return formatter.format(lastActivity);
	}

	//	@OneToMany(mappedBy = "batchUpload", cascade = CascadeType.ALL, fetch = FetchType.LAZY)
	//	@XmlTransient
	//	public List<BatchUploadRow> rows;

	public enum BatchUploadState
	{
		init, remap_columns, remap_columns_submit, remap_entities, remap_entities_submit, preview_browser, preview_browser_submit, finished, error 
	}

	public void refreshLastActivity()
	{
		lastActivity = new Timestamp((new Date()).getTime());
	}

	@SuppressWarnings("unchecked")
	@XmlElement(name="numRows")
	public Long getNumRows()
	{
		if (!Globals.getMarshallingOption(MarshallingOption.EXTENDED_BATCHUPLOAD_INFO))
			return null;
		List<Long> l = Globals.session().createCriteria(BatchUploadRow.class).add(Restrictions.eq("batchUpload", this)).setProjection(Projections.count("id")).list();
		return l.get(0);
	}

	@XmlElement(name="aggregatedRows")
	@SuppressWarnings("unchecked")
	public List<AggregatedRow> getAggregatedRows()
	{
		if (!Globals.getMarshallingOption(MarshallingOption.EXTENDED_BATCHUPLOAD_INFO))
			return null;
		List<AggregatedRow> rows = new ArrayList<AggregatedRow>();
		Criteria c = Globals.session().createCriteria(BatchUploadRow.class)
				.add(Restrictions.eq("batchUpload", this))
				.setProjection(
						Projections.projectionList()
						.add(Projections.groupProperty("type"))
						.add(Projections.groupProperty("status"))
						.add(Projections.count("id")))
						.addOrder(Order.asc("status"))
						.addOrder(Order.asc("type"));
		List<Object[]> l = c.list();
		for (Object[] row : l) 
			rows.add(new AggregatedRow((RowType)row[0], (RowStatus)row[1], (Long)row[2]));

		return rows;
	}

	@XmlElement(name="numRecords")
	@SuppressWarnings("unchecked")
	public Long getNumRecords()
	{
		if (!Globals.getMarshallingOption(MarshallingOption.EXTENDED_BATCHUPLOAD_INFO))
			return null;
		List<Long> l = Globals.session().createCriteria(BatchUploadRow.class)
				.createAlias("ep", "ep")
				.add(Restrictions.eq("batchUpload", this))
				.add(Restrictions.isNull("ep.deleted"))
				.setProjection(Projections.countDistinct("ep.id"))
				.list();
		return l.get(0);
	}

}

@XmlRootElement
class AggregatedRow
{
	@XmlAttribute
	RowType type;
	@XmlAttribute
	RowStatus status;
	@XmlAttribute
	Long count;

	public AggregatedRow()
	{

	}

	public AggregatedRow(RowType type, RowStatus status, Long count)
	{
		this.type = type;
		this.status = status;
		this.count = count;
	}
}