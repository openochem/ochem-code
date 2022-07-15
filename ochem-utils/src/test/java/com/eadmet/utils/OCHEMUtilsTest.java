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

import org.apache.commons.lang.StringUtils;
import org.junit.Assert;
import org.junit.Test;

public class OCHEMUtilsTest
{
	private static final String String = "OCHEM";
	private static final String[] stringArray = new String[]{"OCHEM", "eADMET", "Chemo"};

	@Test
	public void testMySqlCompatibleCompressUncompress() throws IOException
	{
		byte[] compressed = OCHEMUtils.MySqlCompatibleCompress(String);
		byte[] uncompressed = OCHEMUtils.MySqlCompatibleUncompress(compressed);

		Assert.assertArrayEquals(String.getBytes(), uncompressed);
	}

	@Test
	public void testGetCrc32() throws IOException
	{
		long crc32 = OCHEMUtils.getCrc32(String);
		Assert.assertEquals(1858770243l, crc32);
	}

	@Test
	public void testGetMD5() throws IOException
	{
		String md5 = OCHEMUtils.getMD5(String);
		Assert.assertEquals("e82c92a688fa7cf50cdee419f2da41fe", md5);
	}

	@Test
	public void testImplode1() throws IOException
	{
		String imploded = StringUtils.join(Arrays.asList(stringArray), "Informatics");
		Assert.assertEquals("OCHEMInformaticseADMETInformaticsChemo", imploded);
	}

	@Test
	public void testImplode2() throws IOException
	{
		String imploded = StringUtils.join(stringArray, "Informatics");
		Assert.assertEquals("OCHEMInformaticseADMETInformaticsChemo", imploded);
	}
}
