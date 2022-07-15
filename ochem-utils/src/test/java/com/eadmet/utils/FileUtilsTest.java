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

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.junit.Test;

import org.junit.Assert;

public class FileUtilsTest
{
	@Test
	public void fileReadWriteStringTest() throws IOException
	{
		String ps = System.getProperty("file.separator");
		File tmp = FileUtils.createTempDir();
		String referenceString = "ABCDEFGHIJKLMNO";
		String fileName = tmp.getAbsolutePath() + ps + "test.tmp";
		FileUtils.saveStringToFile(referenceString, fileName);
		Assert.assertEquals(referenceString.trim(), FileUtils.getFileAsString(fileName).trim());
	}
	
	@Test
	public void fileReadWriteBytesTest() throws IOException
	{
		String ps = System.getProperty("file.separator");
		File tmp = FileUtils.createTempDir();
		byte[] referenceString = "ABCDEFGHIJKLMNO".getBytes();
		String fileName = tmp.getAbsolutePath() + ps + "test.tmp";
		FileUtils.saveBytesToFile(referenceString, fileName);
		Assert.assertTrue(Arrays.equals(referenceString, FileUtils.getFileAsBytes(fileName)));
	}
	
	@Test
	public void testZip() throws IOException
	{
		File tmp = FileUtils.createTempDir();
		String fileName = tmp.getAbsolutePath() + File.separator + "test.zip";
		File zipFile = new File(fileName);
		
		Assert.assertFalse(zipFile.exists());
		FileUtils.zip(tmp, zipFile);
		
		Assert.assertTrue(zipFile.exists());
		zipFile.delete();
	}
	
}
