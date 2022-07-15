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

import java.io.IOException;
import java.util.Arrays;

import org.junit.Test;

import org.junit.Assert;

public class Base64Test
{
	@Test
	public void testBytes() throws IOException
	{
		byte[] reference = new byte[]{0, 1, 2, 3, 4, 8, 'A'};
		byte[] decoded = Base64.decode(Base64.encodeBytes(reference));
		Assert.assertTrue(Arrays.equals(reference, decoded));
	}
	
	@Test
	public void testBytesApache() throws IOException
	{
		byte[] reference = new byte[]{0, 1, 2, 3, 4, 8, 'A'};
		byte[] decoded = org.apache.commons.codec.binary.Base64.decodeBase64(org.apache.commons.codec.binary.Base64.encodeBase64(reference));
		Assert.assertTrue(Arrays.equals(reference, decoded));
	}
	
	@Test
	public void testCompareBytes() throws IOException
	{
		byte[] reference = new byte[]{0, 1, 2, 3, 4, 8, 'A'};
		byte[] decodedApache = org.apache.commons.codec.binary.Base64.decodeBase64(org.apache.commons.codec.binary.Base64.encodeBase64(reference));
		byte[] decoded = Base64.decode(Base64.encodeBytes(reference));
		
		Assert.assertTrue(Arrays.equals(decodedApache, decoded));	
		
	}
	
}
