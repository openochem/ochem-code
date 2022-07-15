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
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;


@Entity
public class Invite
{
	@Id
	@GeneratedValue
	@Column(name = "invite_id")
	@XmlAttribute
	public Long id;

	@Column
	@XmlTransient
	public String code;

	@XmlTransient
	@ManyToOne
	@JoinColumn(name = "owner_id", referencedColumnName="user_id")
	public User owner;

	@XmlTransient
	@ManyToOne
	@JoinColumn(name = "user_id", referencedColumnName="user_id")
	public User user;

	public static Invite getFree(String name)
	{
		Criteria inviteCriteria = Globals.session().createCriteria(Invite.class)
				.add(Restrictions.eq("code", name))
				.add(Restrictions.isNull("user"));
		@SuppressWarnings("unchecked")
		List<Invite> invites = inviteCriteria.list();
		if (invites.size() > 0)
			return invites.get(0);
		else
			return null;
	}

}
