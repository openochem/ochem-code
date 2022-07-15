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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name= "distribution")
public class ValueDistribution 
{
	public double min;
	public double max;
	
	public double middle;
	
	public List<Double> bins;
	public List<Integer> frequencies;
	
	public void setValues(List<Double> values)
	{
		Collections.sort(values);
		if (values == null || values.size() == 0)
			return;
		min = values.get(0);
		max = values.get(values.size() - 1);
		middle = (min + max) / 2;
		
		bins = new ArrayList<Double>();
		double binWidth = 2.0/3.0 * (max - min) * Math.pow(values.size(), -0.333);
		double bin = min + binWidth;
		while (bin < max)
		{
			bins.add(bin);
			bin += binWidth;
		}
		bins.add(max);
		int[] freq = new int[bins.size()];
		
		for (int i = 0, b = 0; i < values.size(); i++)
		{
			while (bins.get(b) < values.get(i))
				b++;
			freq[b]++;
		}
		
		frequencies = new ArrayList<Integer>();
		for (int i : freq)
			frequencies.add(i);
	}
}
