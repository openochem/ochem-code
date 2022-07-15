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

package qspr.fragmententities;

import java.util.List;

import javax.xml.bind.annotation.XmlType;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
@XmlType(name="FragMapping2" )
public class Mapping2 
{
	public Integer id;

	public String inchi2;

	public String data;

	public Mapping1 mapping1;

	public static Mapping2 get(String fullInchi, String molData)
	{
		return get(fullInchi, molData, false);
	}

	public static Mapping2 get(String inchi1, String inchi2, String molData)
	{
		return get(inchi1, inchi2, molData, false);
	}

	public static Mapping2 get(String fullInchi, String molData, boolean nonFragmentable)
	{
		String[] pieces = fullInchi.split("-");
		if (pieces.length > 1)
			return get(pieces[0], pieces[1], molData, nonFragmentable);
		else
			return get(pieces[0], pieces[0], molData, nonFragmentable);
	}

	@SuppressWarnings("unchecked")
	public static Mapping2 get(String inchi1, String inchi2, String molData, boolean nonFragmentable)
	{
		Mapping2 result = null;
		String cInchi = null;
		if (inchi1.length() == 32) 
		{
			cInchi = inchi1;
		}
		else
		{
			cInchi = inchi1 + "-" + inchi2;
		}

		Criteria criteria = Globals.alternateSession().createCriteria(Mapping2.class)
				.add(Restrictions.eq("inchi2", cInchi));
		List<Mapping2> mapping2s = criteria.list();

		if (mapping2s.size() == 0) // introduce to all
		{

			Mapping1 m1 = (nonFragmentable) ? Mapping1.getNonFragmentable(inchi1, molData) : Mapping1.get(inchi1, molData);

			result = new Mapping2();
			result.mapping1 = m1;
			result.inchi2 = cInchi;
			result.data = molData;
			Globals.alternateSession().saveOrUpdate(result);

		}
		else
			result = mapping2s.get(0);
		return result;
	}

}
