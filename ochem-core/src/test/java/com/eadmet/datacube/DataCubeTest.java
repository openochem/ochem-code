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

package com.eadmet.datacube;

import org.junit.Assert;
import org.junit.Test;

public class DataCubeTest {
	
	@Test
	public void groupingsTest() {
		
		Subkonto sInteger = new Subkonto(Integer.class, "Integer");
		Subkonto sString = new Subkonto(String.class, "String");
		
		DataCube<Double> cube = new DataCube<Double>(Double.class);
		cube.setupGroupingOrder(sInteger, sString);
		cube.addValue(1.0, "AAA", 1);
		cube.addValue(2.0,  "BBB", 1);
		cube.addValue(4.0,  2, "BBB");

		Assert.assertEquals(7.0, cube.getValue(), 0.01);
		Assert.assertEquals(3.0, cube.getValue(1), 0.01);
		Assert.assertEquals(4.0, cube.getValue(2), 0.01);
	}
}