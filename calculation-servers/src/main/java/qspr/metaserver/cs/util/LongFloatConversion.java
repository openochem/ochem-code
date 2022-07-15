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

package qspr.metaserver.cs.util;

import java.io.IOException;

public class LongFloatConversion {

	static long PRECISION = 1000;

	public static long [] getInt(long l) throws IOException {
		long vals[] = new long[3]; 

		if(l>PRECISION*PRECISION*PRECISION)throw new IOException("Provided id is too large");
		vals[0] = l % PRECISION;
		vals[1] = (l/PRECISION)%PRECISION;
		vals[2] = (l/(PRECISION*PRECISION))%PRECISION;

		check(vals, l);

		return vals;
	}

	public static long getLong(long vals[]) {
		long s = 0, n =1;
		for(int i=0;i<vals.length;i++) {
			s += vals[i]*n;
			n *= PRECISION;
		}
		return s;
	}

	static void check(long vals[],long l) throws IOException {
		if(l != getLong(vals)) throw new IOException("Values are different ones: " + l + "  " + getLong(vals));
	}

}
