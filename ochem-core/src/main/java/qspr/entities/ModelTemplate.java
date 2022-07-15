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

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.metaserver.protocol.Task;
import qspr.workflow.utils.QSPRConstants;

@Entity
@XmlRootElement(name = "model-template")
public class ModelTemplate 
{
	@Id
	@Column(name = "model_template_id")
	@GeneratedValue(strategy = GenerationType.AUTO)
	@XmlAttribute
	Long id;

	@Column
	@XmlAttribute
	public String name;

	@Column(name = "displayed_name")
	@XmlElement
	public String displayedName;

	@ManyToOne(cascade={CascadeType.ALL})
	@JoinColumn(name = "teacher_task_template")
	@XmlTransient
	public Attachment<Task> teacherTaskTemplate;

	@ManyToOne(cascade={CascadeType.ALL})
	@JoinColumn(name = "applier_task_template")
	@XmlTransient
	public Attachment<Task> applierTaskTemplate;

	@Column(name = "is_support_multilearning")
	@XmlAttribute
	public Boolean is_Support_Multilearning;

	@Column(name = "is_suggested_model")
	@XmlAttribute
	public Boolean isSuggestedModel;

	@Column
	public boolean invisible;

	public ModelTemplate()
	{

	}

	public ModelTemplate(String name)
	{
		this.name = name;
	}

	@XmlAttribute
	public Boolean getUploadable()
	{
		if (QSPRConstants.MLRA.equals(name))
			return true;
		return null;
	}

	public boolean isDescriptorCalculationOnly()
	{
		return name.equals("Descriptors Without Model");
	}

}
