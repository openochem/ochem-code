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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.hibernate.Criteria;
import org.hibernate.criterion.Order;

public class CriteriaWrapper
{
	public Criteria criteria;
	public List<Order> order = new ArrayList<Order>();

	public Map<String, String> aliases = new LinkedHashMap<String, String>();

	public CriteriaWrapper(Criteria criteria)
	{
		this.criteria = criteria;
	}

	public Criteria createAlias(String path, String alias)
	{
		if (!aliases.containsKey(alias))
		{
			aliases.put(alias, path);
			criteria.createAlias(path, alias);
		}
		return criteria;
	}
	
	public void addAliasesTo(Criteria c) {
		for (String alias : aliases.keySet())
			c.createAlias(aliases.get(alias), alias);
	}
	
	public void addOrderTo(Criteria c) {
		for (Order o : order)
			c.addOrder(o);
	}

	public boolean hasAlias(String alias)
	{
		return aliases.containsKey(alias);
	}
}