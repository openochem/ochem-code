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
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;

import qspr.Globals;

@Entity
public class UserAction {
	
	@Id
	@GeneratedValue
	@Column(name = "useraction_id")
	public Long id;

	@ManyToOne
	@JoinColumn(name = "action_id")
	public Action action;

	@ManyToOne
	@JoinColumn(name = "user_id")
	public User user;
	
	public static void save(User user, Action action)
	{
		UserAction ua = new UserAction();
		ua.action = action;
		ua.user = user;
		Globals.session().save(ua);	
	}
}
