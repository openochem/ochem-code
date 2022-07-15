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

import org.junit.Test;

import org.junit.Assert;

public class CryptUtilsTest
{
	public static final String REFERENCE_STRING = "ABCDEFG";
	@Test
	public void encodeDecodeTest() throws Exception
	{
		
		Assert.assertEquals(REFERENCE_STRING, CryptUtils.desDecode(CryptUtils.desEncode(REFERENCE_STRING)));
	}
	
	@Test
	public void testConsistency() throws Exception
	{
		Assert.assertEquals("eADkA4rb/XA=", CryptUtils.desEncode(REFERENCE_STRING));
	}
}
