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

package qspr.schedule;

import qspr.tests.PeriodicTestRunner;

import com.eadmet.utils.config.ConfigurableClass;
import com.eadmet.utils.config.ConfigurableProperty;

// Just a wrapper for the Quartz scheduler
@GlobalMaintenanceJob
@ConfigurableClass(name = "ochem")
public class PeriodicTestsTask extends OchemCronjobTask
{
	@ConfigurableProperty(name = "periodic_tests", comment = "When enabled, the server will periodically run all the tests (hourly, bihourly or daily, depending on the test type)")
	static public boolean runPeriodicTests = true;
	
	@ConfigurableProperty(name = "periodic_tests_type", comment = "Which periodic tests should run?")
	static public String periodicTestsType;
	
	@Override
	public boolean shouldRun()
	{
		return super.shouldRun() && PeriodicTestsTask.runPeriodicTests;
	}
	
	@Override
	public void executeTask() throws Exception
	{
		PeriodicTestRunner runner = new PeriodicTestRunner();
		runner.filter.runSpecificTestType = periodicTestsType;
		runner.run();
	}
}

