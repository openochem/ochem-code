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

import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

import qspr.metaserver.util.aggregator.AggregatedValue;

public class AggregatedValueTest
{
	@Test
	public void basicTest()
	{
		Random r = new Random();
		AggregatedValue av = new AggregatedValue(true);
		for (int i = 0; i < 10000; i++)
			av.addValue(r.nextGaussian());
		
		Assert.assertEquals(0.0, av.getMean(), 0.05);
		Assert.assertEquals(1.96, av.get95ConfidenceInterval(), 0.07);
		Assert.assertEquals(1.0, av.getSTD(), 0.05);
	}
}
