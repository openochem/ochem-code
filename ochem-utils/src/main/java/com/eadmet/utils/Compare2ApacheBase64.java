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

import java.util.Arrays;

import org.apache.commons.codec.binary.Base64;

import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;


public class Compare2ApacheBase64
{

	public static final String EMAIL_ADMIN = "itetko@vcclab.org";
	
	public static void compareEncode(byte[] cipherText, String encodedCipherText)
	{
		String apache = Base64.encodeBase64String(cipherText);
		if ( ! apache.equals(encodedCipherText))
			Mailer.postMailAsync(new Email(EMAIL_ADMIN, "Base64 comparison failed", "Comparison of base64 encoding failed.\n "));
	}
	
	public static void compareDecode(String arg, byte[] decodedCipher)
	{
		byte[] apache = Base64.decodeBase64(arg);
		if ( ! Arrays.equals(decodedCipher, apache))
			Mailer.postMailAsync(new Email(EMAIL_ADMIN, "Base64 comparison failed", "Comparison of base64 decoding failed.\n "));
	}

	
	
}
