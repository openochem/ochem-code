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
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;

import qspr.workflow.utils.QSPRConstants;

@Entity
public class ArticleUserPdf 
{
	static public final int PDF = 1;
	static public final int EXCEL = 2;
	static public final int SDF = 3;
	
	@Id
	@GeneratedValue
	@Column(name = "articleuserpdf_id")
	public Integer id;
	
	@ManyToOne
	@JoinColumn(name = "article_id")
	public Article article;
	
	@ManyToOne
	@JoinColumn(name = "user_id")
	public User user;
	
	@ManyToOne(fetch = FetchType.LAZY, optional = false)
	@JoinColumn(name = "attachment_id")
	public Attachment attachedFile;
	
	@Column
	@XmlAttribute
	public Integer type;
	
	public String toString()
	{
		return type == ArticleUserPdf.PDF ? QSPRConstants.PDF : QSPRConstants.EXCEL;
	}
	
}
