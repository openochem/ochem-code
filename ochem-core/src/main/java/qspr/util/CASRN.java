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

package qspr.util;

public class CASRN
{
	public static String checkCasrnSyntax(String s)
	{
		
		try {
			// 												   containing any word-character										  
			if (s == null || s.equals("") || s.equals("0") || s.equals("*") || s.matches(".*[a-zA-Z].*"))
				return null;
			int i;
			for (i = 0; s.charAt(i) == '0'; i++);
			if (i > 0)
				s = s.substring(i);
			if (s.length() == 0)
				return null;
			int j, n1 = s.indexOf("-"), n = s.length(), n2 = s.lastIndexOf("-"), sum;
			if (n1 != -1 && (n1 != n - 5 || n2 != n - 2))
				return null;
			for (i = n1 == -1 ? n - 2 : n - 3, j = 1, sum = 0; i >= 0; i--)
				if (i != n1)
					sum += Integer.parseInt(s.substring(i, i + 1)) * j++;
			j = Integer.parseInt(s.substring(n - 1));
			if ((sum - j) % 10 != 0)
				return null;
			if (n1 == -1)
				s = s.substring(0, n - 3) + "-" + s.substring(n - 3, n - 1) + "-" + j;
			return s;
		} catch (NumberFormatException e) {
			// TODO: handle exception
		} catch (Exception ea) {
			return null;
		}
		return null;
	}
	
}
