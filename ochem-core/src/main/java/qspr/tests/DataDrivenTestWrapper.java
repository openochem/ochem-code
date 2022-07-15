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

package qspr.tests;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.entities.Session;
import qspr.entities.User;
import qspr.util.WrapperThread;

public abstract class DataDrivenTestWrapper extends WrapperThread {

	public abstract void wrappedTest() throws Exception;
	protected String prefix;
	public boolean cleanupBefore = true;
	
	public DataDrivenTestWrapper(String prefix)
	{
		this.prefix = prefix;
	}
	
	public DataDrivenTestWrapper setCleanupBefore(boolean cleanupBefore)
	{
		this.cleanupBefore = cleanupBefore;
		return this;
	}
	
	@Override
	public void wrapped() throws Exception 
	{ 
		if (cleanupBefore)
			OCHEMTestHelpers.cleanupByPrefix(prefix);
		
		ThreadScope.get().disableTrackChanges = true;
		Session userSession = OCHEMTestHelpers.generateRandomTestUser(prefix);
		
		Globals.session().flush();
		ThreadScope.get().userSession = userSession;
		Globals.restartAllTransactions(true);
		try
		{
			wrappedTest();
		} finally
		{
			// make sure you have a active transaction to clean up
			Globals.commitAllTransactions();
			Globals.startAllTransactions();
			
			User user = User.getByLogin(userSession.user.login);
			OCHEMTestHelpers.cleanup(user, true);
			Globals.restartAllTransactions(true); //Individual cleanup
			OCHEMTestHelpers.cleanupByPrefix(prefix);
			Globals.restartAllTransactions(true); //Potential former crashes cleanup
		}
	}
	
}
