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

import junit.framework.JUnit4TestAdapter;
import junit.framework.TestSuite;

import org.junit.runner.RunWith;
import org.junit.runner.manipulation.NoTestsRemainException;
import org.junit.runners.AllTests;

import qspr.Globals;

/**
 * Allows to run all the periodic tests as a jUnit test from Eclipse
 * @author midnighter
 *
 */
@RunWith(AllTests.class)
public final class PeriodicTestsSuite {

	@SuppressWarnings("rawtypes")
	public static TestSuite suite() {
		TestSuite suite = new TestSuite();
		//Globals.executableDirectory="/ews/ochem-webapp/src/main/webapp/WEB-INF/executables";
		Globals.commonDownloadDirectory = "/ews/ochem-webapp/src/main/webapp/WEB-INF/exports";

		PeriodicTestsListener listener = new PeriodicTestsListener(null, true);
		PeriodicTestFilter filter = new PeriodicTestFilter(listener);
		filter.disableScheduler = true;
		//filter.runSpecificTestType = "general";
		listener.saveToDB = false;
		for (Class clazz : PeriodicTestRunner.annotatedTestClasses)
		{
			JUnit4TestAdapter adapter = new JUnit4TestAdapter(clazz);
			try
			{
				adapter.filter(filter);
			} catch (NoTestsRemainException e)
			{
			}
			suite.addTest(adapter);
		}

		return suite;
	}
}
