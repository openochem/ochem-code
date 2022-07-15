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

package qspr;

import org.junit.Assert;
import org.junit.Test;

public class OCHEMConfigurationTest
{
	@Test
	public void ipAddressTest() {
		OCHEMConfiguration.trustedIPAddressesString = "123.456, 123.666";
		OCHEMConfiguration.addTrustedAddress("888");
		Assert.assertTrue(OCHEMConfiguration.isTrustedAddress("123.456.777.777"));
		Assert.assertTrue(OCHEMConfiguration.isTrustedAddress("123.888.777.777"));
		Assert.assertFalse(OCHEMConfiguration.isTrustedAddress("123.889.777.777"));
	}
}
