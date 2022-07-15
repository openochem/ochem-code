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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "short-condition")
public class ShortCondition implements Serializable
{
	private static final long serialVersionUID = 1L;

	public long id; // This is currently used for condition_id
	public Double defaultValue = 0.;
	public long unitId;

	// Mergings of options of qualitative conditions
	public List<OptionsMerging> optionsMergings;

	// an aid structure
	protected HashMap<Long, Long> mergingEndpoints;

	public ShortCondition(ShortCondition desc) {
		id = desc.id;
		defaultValue = desc.defaultValue;
		unitId = desc.unitId;
		optionsMergings = desc.optionsMergings;
	}

	
	public ShortCondition()
	{
		// To make JAXB happy
	}

	public ShortCondition(long id)
	{
		this.id = id;
	}

	public void addMerging(long id1, long id2)
	{
		// Do not distinguish the two specified options of this condition
		if (optionsMergings == null)
			optionsMergings = new ArrayList<OptionsMerging>();
		OptionsMerging om = new OptionsMerging();
		om.id1 = id1;
		om.id2 = id2;
		optionsMergings.add(om);
	}

	public long getAnalogueOption(long id)
	{
		Long res = mergingEndpoints.get(id);
		return res == null ? id : res;
	}

	public void calculateMergings()
	{
		mergingEndpoints = new HashMap<Long, Long>();

		if (optionsMergings == null)
			return;

		for (OptionsMerging om : optionsMergings) {
			mergingEndpoints.put(om.id1, om.min());
			mergingEndpoints.put(om.id2, om.min());
		}

		boolean more = true;
		while (more)
		{
			more = false;

			// For every connection
			for (OptionsMerging om : optionsMergings) 
			{
				long curMp;
				long mp1 = mergingEndpoints.get(curMp = mergingEndpoints.get(om.id1));
				long mp2 = mergingEndpoints.get(mergingEndpoints.get(om.id2));
				if (mp1 != mp2 || curMp != mp1)
				{
					long min = Math.min(mp1, mp2);
					mergingEndpoints.put(om.id1, min);
					mergingEndpoints.put(om.id2, min);
					more = true;
				}
			}
		}
	}


	// Merging of two options of a qualitative condition
	public static class OptionsMerging implements Serializable
	{
		private static final long serialVersionUID = 1L;

		public long id1;
		public long id2;

		public long min()
		{
			return Math.min(id1, id2);
		}
	}

}
