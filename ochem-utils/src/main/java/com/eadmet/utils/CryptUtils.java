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

import java.security.MessageDigest;
import java.util.Arrays;

import javax.crypto.Cipher;
import javax.crypto.SecretKey;
import javax.crypto.spec.IvParameterSpec;
import javax.crypto.spec.SecretKeySpec;

public class CryptUtils 
{
	private final static String masterPassword = "masterPasswordOCHEM";

	public static String desEncode(String arg) throws Exception
	{
		MessageDigest md = MessageDigest.getInstance("md5");
		byte[] digestOfPassword = md.digest(masterPassword.getBytes("utf-8"));
		byte[] keyBytes = Arrays.copyOf(digestOfPassword, 24);
		for (int j = 0,  k = 16; j < 8;)
			keyBytes[k++] = keyBytes[j++];

		SecretKey key = new SecretKeySpec(keyBytes, "DESede");
		IvParameterSpec iv = new IvParameterSpec(new byte[8]);
		Cipher cipher = Cipher.getInstance("DESede/CBC/PKCS5Padding");
		cipher.init(Cipher.ENCRYPT_MODE, key, iv);

		byte[] plainTextBytes = arg.getBytes("utf-8");
		byte[] cipherText = cipher.doFinal(plainTextBytes);
		String encodedCipherText = Base64.encodeBytes(cipherText);
		Compare2ApacheBase64.compareEncode(cipherText, encodedCipherText);
		return encodedCipherText;
	}

	public static String desDecode(String arg) throws Exception
	{
		MessageDigest md = MessageDigest.getInstance("md5");
		byte[] digestOfPassword = md.digest(masterPassword.getBytes("utf-8"));
		byte[] keyBytes = Arrays.copyOf(digestOfPassword, 24);
		for (int j = 0,  k = 16; j < 8;)
			keyBytes[k++] = keyBytes[j++];

		SecretKey key = new SecretKeySpec(keyBytes, "DESede");
		IvParameterSpec iv = new IvParameterSpec(new byte[8]);
		Cipher decipher = Cipher.getInstance("DESede/CBC/PKCS5Padding");
		decipher.init(Cipher.DECRYPT_MODE, key, iv);
		byte[] decodedCipher = Base64.decode(arg);
		Compare2ApacheBase64.compareDecode(arg, decodedCipher);
		byte[] plainText = decipher.doFinal(decodedCipher);
		return new String(plainText);     
	}

	public static boolean isMasterPassword(String password)
	{
		return password.equals(masterPassword);
	}

	public static void main(String... args) throws Exception
	{
		System.out.println(desEncode("demodata"));
		System.out.println(desDecode("Hp6GCS0TaVworus2nIY/dw=="));
	}
}
