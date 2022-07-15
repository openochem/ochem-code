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

package qspr.dao;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.Globals;
import qspr.entities.Session;
import qspr.entities.User;

public class UserDAOImpl implements UserDAO{

	@Override
	public User getById(long id) {
		return (User) Globals.session().get(User.getCurrentClass(), id);
	}

	/**
	 * 
	 * @param rows
	 * @param multiplier
	 */

	@Override
	public void checkEligibility(long rows, int multiplier)throws UserFriendlyException{

	}

	@Override
	public User getSessionUser() {
		return Globals.userSession().user;
	}

	@Override
	public User getBySessionId(long id) {
		Session session = (Session) Globals.session().get(Session.class, id);
		return getById(session.user.id);
	}

}
