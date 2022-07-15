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

package qspr.entities;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

@XmlRootElement(name = "model-microattachment")
public class ModelMicroAttachment 
{
	@XmlTransient
	public Set<Integer> excludedBasketEntries = new HashSet<Integer>();
	
	@XmlElement(name="excluded-basketentries")
	public String getExcludedBasketEntries()
	{
		if (excludedBasketEntries.size() == 0)
			return "";
		
		List<Integer> ids = new ArrayList<Integer>();
		ids.addAll(excludedBasketEntries);
		String s = ids.get(0).toString();
		
		for (int i=1; i<ids.size(); i++)
			s+=(","+ids.get(i));
		
		return s;
	}
	
	public void setExcludedBasketEntries(String s)
	{
		excludedBasketEntries.clear();
		String[] entries = s.split(",");
		for (String entry : entries)
			if (!entry.equals(""))
				excludedBasketEntries.add(Integer.valueOf(entry));
	}

}
