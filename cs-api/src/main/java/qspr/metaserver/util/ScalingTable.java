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

package qspr.metaserver.util;

import java.io.IOException;

// To clarify: is there any class implementing this interface except of PrototypeTable?
// if not, these functions can be moved to PrototypeTable
public interface ScalingTable {

	public int getColumnSize();
	
	public double[] getMin(int size) throws IOException;
	public double[] getMax(int size) throws IOException;
	public double[] getStd(int size) throws IOException;
	public double[] getMean(int size) throws IOException;

}
