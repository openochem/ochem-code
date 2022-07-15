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

package qspr.services;

import java.sql.Timestamp;
import java.util.Calendar;

import qspr.Globals;
import qspr.entities.Session;
import qspr.entities.User;

import com.eadmet.exceptions.UserFriendlyException;

/**
 * A set of simple services to work with the model
 * @author itetko
 *
 */

public class LoginService {

	public boolean useApplierCache = true;

	public String login(String username, String password) {
		try {
			Globals.startAllTransactions();
			
			User user = null;
			if (!"guest".equals(username))
			{
				
				user = User.getByLogin(username);
	
				if (user == null)
					throw new UserFriendlyException("Unknown user");
	
				if (!user.authorizeWithSecret(password))
					throw new UserFriendlyException("Wrong password");
			}

			Session newUserSession = Session.createNewSession(user);
			newUserSession.ipAddress = "WebService";
			newUserSession.time = new Timestamp(Calendar.getInstance().getTimeInMillis());
			Globals.session().save(newUserSession);

			Globals.commitAllTransactions(); // FIXME: Do not start a fragments-DB transaction

			return newUserSession.guid;
		} catch (UserFriendlyException e) {
			Globals.rollbackAllTransactions();
			throw e;
		} catch (Exception e) {
			e.printStackTrace();
			throw new UserFriendlyException(e.getMessage());
		}
	}
	
}
