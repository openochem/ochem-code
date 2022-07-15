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

import java.util.ArrayList;
import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.ThreadScope;

@Entity(name = "Gruppe")
@XmlRootElement(name = "group")
public class Group 
{

	@Id
	@Column(name = "group_id")
	@GeneratedValue
	@XmlAttribute
	public Long id;

	@Column
	@XmlElement(name = "name")
	public String name;

	@XmlElement(name = "member")
	public List<User> getUsers()
	{
		if (ThreadScope.get().controller.equals("epbrowser"))
		{

			@SuppressWarnings("unchecked")
			List<User> users = Globals.session().createCriteria(User.getCurrentClass()).add(Restrictions.eq("group", this)).list();
			List<User> result = new ArrayList<User>();
			for (User user : users) {
				// Prevent JAXB looping this way
				User usr = User.getNewInstance();
				usr.id = user.id;
				usr.login = user.login;
				result.add(usr);
			}

			return result;
		}
		return null;
	}


	public String toString()
	{
		return name;
	}

	public boolean equals(Object obj)
	{
		Group group = (Group) obj;
		return (group == null) ? false: ((id == null) ? (group.id == null) : id.equals(group.id));
	}

}
