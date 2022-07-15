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
import cern.colt.Timer;

public class TestTheTest {
	public static void main(String[] args) {
		//Globals.executableDirectory="/ews/ochem-webapp/src/main/webapp/WEB-INF/executables";
		Globals.commonDownloadDirectory = "/ews/ochem-webapp/src/main/webapp/WEB-INF/exports";
		PeriodicTestRunner runner = new PeriodicTestRunner();

//		runner.filter.runSpecificTestType = "selenium";
//		runner.filter.runSpecificTest = "BasketBrowser"; 
//		runner.filter.runSpecificTest = "modelConsistency";
		
		runner.filter.runSpecificTest = "newBatchUploadTest"; 

		runner.filter.disableScheduler = true;
		runner.listener.saveToDB = false;
		//runner.saveReport = true;
		Timer t = new Timer(); t.start();
		runner.run();
		System.out.println("Tests took: " + t.seconds() + "secs.");
	}
}

