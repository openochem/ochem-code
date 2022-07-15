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

public class Mathematics {

	// return phi(x) = standard Gaussian pdf
	public static double phi(double x) {
	    return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
	}
	
	// return phi(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma
	public static double phi(double x, double mu, double sigma) {
	    return phi((x - mu) / sigma) / sigma;
	}
	
	// return Phi(z) = standard Gaussian cdf using Taylor approximation
	public static double Phi(double z) {
	    if (z < -8.0) return 0.0;
	    if (z >  8.0) return 1.0;
	    double sum = 0.0, term = z;
	    for (int i = 3; sum + term != sum; i += 2) {
	        sum  = sum + term;
	        term = term * z * z / i;
	    }
	    return 0.5 + sum * phi(z);
	}
	
	// return Phi(z, mu, sigma) = Gaussian cdf with mean mu and stddev sigma
	public static double Phi(double z, double mu, double sigma) {
	    return Phi((z - mu) / sigma);
	} 
	
	// Compute z such that Phi(z) = y via bisection search
	public static double PhiInverse(double y) {
	    return PhiInverse(y, .00000001, -8, 8);
	} 
	
	// bisection search
	private static double PhiInverse(double y, double delta, double lo, double hi) {
	    double mid = lo + (hi - lo) / 2;
	    if (hi - lo < delta) return mid;
	    if (Phi(mid) > y) return PhiInverse(y, delta, lo, mid);
	    else              return PhiInverse(y, delta, mid, hi);
	}

}
