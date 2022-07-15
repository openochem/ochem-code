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

import java.util.UUID;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;

/**
 * A "model identity" is an alias for a model distributed to external parties. 
 * It allows to :
 *  (a) transparently substitute a model without notifying the external party.
 *  (b) track model usage (e.g., via web services) by different parties.
 *  
 *  As discussed on the development meeting on 3rd June 2013.
 *  Originally inspired by the potential collaboration with Chemicalize
 * 
 * @author midnighter
 *
 */
@Entity
@XmlRootElement(name = "model-identity")
public class ModelIdentity 
{
	@Id
	@GeneratedValue
	@Column(name = "mi_id")
	@XmlAttribute
	protected Integer id;
	
	/**
	 * This GUID is the identity itself.
	 * External parties will refer to models via GUIDs
	 */
	@Column
	@XmlAttribute
	public String guid;
	
	@Column
	@XmlAttribute
	public String comment;

	/**
	 * The actual model hidden behind this identity
	 */
	@ManyToOne
	@JoinColumn(name = "model_id")
	public Model model;
	
	/**
	 * Get an identity by its GUID
	 */
	public static ModelIdentity getByGUID(String guid)
	{
		return (ModelIdentity) Globals.session().createCriteria(ModelIdentity.class).add(Restrictions.eq("guid", guid)).uniqueResult();
	}
	
	public ModelIdentity(Model m)
	{
		this.model = m;
		
		// Generate a new random GUID
		this.guid = UUID.randomUUID().toString();
	}
	
	public ModelIdentity()
	{
		
	}
}
