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

package com.eadmet.utils;

import org.junit.Assert;

import org.junit.Test;

public class NumericalValueStandardizerTest
{
	@Test
	public void basicTest()
	{
		Assert.assertEquals("0.009", NumericalValueStandardizer.getSignificantDigitsStr(0.0088, 1));
		Assert.assertEquals("0.01", NumericalValueStandardizer.getSignificantDigitsStr(0.0099, 1));
		Assert.assertEquals("0.3", NumericalValueStandardizer.getSignificantDigitsStr(1.0/3, 1));
		Assert.assertEquals("0.33", NumericalValueStandardizer.getSignificantDigitsStr(1.0/3, 2));
		Assert.assertEquals("120000", NumericalValueStandardizer.getSignificantDigitsStr(123456.1234, 2));
		Assert.assertEquals("123000", NumericalValueStandardizer.getSignificantDigitsStr(123456.1234, 3));
		Assert.assertEquals("1.2E8", NumericalValueStandardizer.getSignificantDigitsStr(123456789, 2));
		
	}
}
