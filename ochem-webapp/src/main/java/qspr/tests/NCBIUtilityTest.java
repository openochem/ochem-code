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

import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;

import javax.xml.bind.JAXBException;

import org.junit.Test;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.entities.Molecule;
import qspr.util.NCBI_Utility;

@PeriodicTest(schedule = ScheduleType.DAILY, type = "general")
public class NCBIUtilityTest
{

	@Test(timeout = 180000)
	public void testGetMoleculeByName() throws MalformedURLException, UnsupportedEncodingException, JAXBException
	{
		if(OCHEMConfiguration.allowExternalServices){
			Globals.startAllTransactions();
			Molecule molecule = NCBI_Utility.getMoleculeByName("aspirine");
			assertNotNull(molecule);
			Globals.rollbackAllTransactions();
		}
	}

	@Test(timeout = 180000)
	public void testGetMoleculeByNameFailure() throws MalformedURLException, UnsupportedEncodingException, JAXBException
	{
		if(OCHEMConfiguration.allowExternalServices){
			Globals.startAllTransactions();
			Molecule molecule = NCBI_Utility.getMoleculeByName("Strandkorb");
			assertNull(molecule);
			Globals.rollbackAllTransactions();
		}
	}

}
