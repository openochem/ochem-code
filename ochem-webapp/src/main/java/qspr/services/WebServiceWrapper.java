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

import qspr.Globals;
import qspr.LoadLogger;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.entities.Session;
import qspr.entities.User;
import qspr.metaserver.protocol.Task;
import qspr.toxicity.DispatcherServletWrapper;
import qspr.util.DynaWrap;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;

//TODO: Refactoring required - most (all?) of the web service methods should be called via this wrapper
public abstract class WebServiceWrapper {
	protected ModelResponse response = new ModelResponse();

	public abstract void run() throws Exception;

	public String sessionGUID;
	public boolean throwExceptions = false;
	protected Session session;

	public ModelResponse execute() {
		try {
			if (OCHEMConfiguration.logLoad)
				LoadLogger.log();
			try {
				if(!Globals.areTransactionsRunning())
					Globals.startAllTransactions();
			}
			catch(java.lang.Throwable e) {
			}

			ThreadScope.get().threadClassLoader = DispatcherServletWrapper.class.getClassLoader();
			ThreadScope.get().updateSessionTime = false;
			if (sessionGUID != null)
				ThreadScope.get().userSession = session = Session.getByGUID(sessionGUID);
			else
				ThreadScope.get().userSession = session = getWebserviceUserSession();

			if (ThreadScope.get().userSession == null)
				ThreadScope.get().userSession = session = getWebserviceUserSession();
			//				throw new UserFriendlyException("Unknown session GUID");
			run();
			Globals.commitAllTransactions();
		} catch (UserFriendlyException e) {
			e.printStackTrace();
			Globals.rollbackAllTransactions();
			response.setStatus(Task.ERROR);
			response.setDetailedStatus("Error: " + e.getMessage());
		} catch (Exception e) {
			e.printStackTrace();
			Globals.rollbackAllTransactions();
			response.setStatus(Task.ERROR);
			response.setDetailedStatus("Error: " + OCHEMUtils.exceptionToString(e));
		} finally {
			ThreadScope.reset();
		}

		return response;
	}

	public WebServiceWrapper() {

	}

	private Session getWebserviceUserSession() throws Exception
	{
		Session session = Session.getFirstSession("test");
		if (session == null)
		{
			User user = User.getByLogin("test");
			if (user == null)
			{
				user = User.getNewInstance();
				user.login = "test";
				if (user.isExtended()) {
					DynaWrap extended = user.dynaWrapped();
					User.getExtended().getMethod("setPassword", String.class).invoke(user, "t3stus3r");
					extended.setField("firstName", "Test");
					extended.setField("organisation", "Academic");
				}
				user.rank = User.NOTVALIDATED;
				user.licenseVersion = Globals.CURRENT_LICENSE_VERSION;
				Globals.session().save(user);
				session = Session.createNewSession(user);
				Globals.session().save(session);
			}
		}

		return session;
	}


	public WebServiceWrapper(String sessionGUID, boolean throwExceptions) 
	{
		this.sessionGUID = sessionGUID;
		this.throwExceptions = throwExceptions;
	}
}
