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

import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;

@Entity
@XmlRootElement(name = "descriptor")
public class Descriptor 
{
	@Id
	@GeneratedValue
	@Column(name = "descriptor_id")
	@XmlAttribute
	public Long id;

	@Column
	@XmlAttribute
	public String name;

	@Column
	@XmlAttribute
	public String program;

	@XmlAttribute
	@Transient
	public Boolean selected;

	@Transient
	@XmlAttribute
	public Double factor;

	@SuppressWarnings("unchecked")
	public static Descriptor getByName(String _name)
	{
		List<Descriptor> condList = Globals.session().createCriteria(Descriptor.class).add(Restrictions.eq("name", _name)).list();
		if (condList.size() > 0)
			return condList.get(0);
		else
			return null;
	}
}
