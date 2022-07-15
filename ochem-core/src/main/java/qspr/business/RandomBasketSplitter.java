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

package qspr.business;
import java.util.ArrayList;
import java.util.List;


public class RandomBasketSplitter extends BasketSplitter
{
	public int validationSetPercentage = 20;
	
	@Override
	protected List<Integer> getSplittedIndices(List<Long> epIDs) 
	{
		List<Integer> selectedIndices = new ArrayList<Integer>();
		for (int i = 0; i < epIDs.size(); i++)
			selectedIndices.add(i);
		
		long selectionSize = Math.round(1.0 * epIDs.size() * validationSetPercentage / 100); 
		for (int i = 0; i < selectionSize; i++)
		{
			int s = Long.valueOf(Math.round((float) 1.0f * Math.random() * selectedIndices.size())).intValue();
			if (s == selectedIndices.size()) // it happens
				s--;
			selectedIndices.remove(s);
		}
		
		return selectedIndices;
	}
}