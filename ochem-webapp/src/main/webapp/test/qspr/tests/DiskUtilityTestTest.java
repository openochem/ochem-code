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

import java.sql.DriverManager;

import javax.mail.MessagingException;

import org.junit.Test;

import qspr.workflow.utils.QSPRConstants
;

import com.eadmet.utils.Mailer;

public class DiskUtilityTestTest {

	@Test
	public void testname() throws Exception {
		
		Mailer.postMail(QSPRConstants.EMAIL_ADMIN, "DiskUtilityTestTest is in use!!!", "Please do not delete me");
		
		TestResult dsTest = DiskUtilityTest.diskSpaceTest();

		try {
			if ( ! dsTest.succeeded) {
				dsTest.postMail(QSPRConstants.EMAIL_ADMIN, "Warning: "+dsTest.name+" has less then 10 GB space left!", dsTest.detailedStatus);
			}	
		} catch (MessagingException e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) throws Exception {
		DiskUtilityTestTest d = new DiskUtilityTestTest();
		d.testname();
	}
	
}
