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

import org.junit.rules.TestRule;
import org.junit.runner.Description;
import org.junit.runners.model.Statement;

import qspr.Globals;
import qspr.metaserver.transport.TransportFactory;

/**
* A test rule that:
* - manages DB transactions
* - creates a test user 
* - runs all necessary cleanups afterwards
* 
*  It is based on the DataDrivenTestWrapper class
*  
* @author robert, midnighter
*/
public class LoginRule implements TestRule
{
	private String testPrefix;
	public boolean keepTransactionOpen = false;
	public boolean cleanupBefore = true;
	
	public LoginRule(String testPrefix, boolean keepTransactionOpen)
	{
		this.testPrefix = testPrefix;
		this.keepTransactionOpen = keepTransactionOpen;
	}
	
	public LoginRule setCleanupBefore(boolean cleanupBefore)
	{
		this.cleanupBefore = cleanupBefore;
		return this;
	}
	
	@Override
	public Statement apply(final Statement statement, Description description)
	{
		return new Statement()
		{
			@Override
			public void evaluate() throws Throwable
			{
				DataDrivenTestWrapper wrapper = new DataDrivenTestWrapper(testPrefix)
				{
					@Override
					public void wrappedTest() throws Exception
					{
						try
						{
							if (!keepTransactionOpen)
								Globals.commitAllTransactions();
							statement.evaluate();
						}
						catch (Exception e)
						{
							throw e;
						}
						catch (Throwable e)
						{
							throw new RuntimeException(e);
						}
						finally
						{
							TransportFactory.clearThreadTransport();
							if (!keepTransactionOpen)
								Globals.startAllTransactions();
						}
					}
				};
				
				wrapper.setCleanupBefore(cleanupBefore);
				
				wrapper.run();
				if (wrapper.exception != null)
					throw wrapper.exception;
			}
		};
	}
	
}