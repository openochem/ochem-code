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

import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;

@Entity
@XmlRootElement(name = "predicate")
public class Predicate 
{

	@Id
	@GeneratedValue
	@Column(name = "predicate_id")
	@XmlAttribute
	public Long id;

	@Column(name = "short_name")
	@XmlAttribute
	public String shortName;

	@Column
	@XmlAttribute
	public String name;

	@SuppressWarnings("unchecked")
	public static Predicate get(String shortName)
	{
		List<Predicate> list = Globals.session().createCriteria(Predicate.class)
				.add(Restrictions.eq("shortName", shortName)).list();
		if (list.size() > 0)
			return list.get(0);

		list = Globals.session().createCriteria(Predicate.class)
				.add(Restrictions.eq("name", shortName)).list();

		if (list.size() > 0)
			return list.get(0);

		return null;
	}

	public boolean isOrdering()
	{
		return shortName.startsWith(">") || shortName.startsWith("<");
	}

	public boolean isInterval()
	{
		return shortName.equals("-");
	}

	public boolean isAccuracy()
	{
		return shortName.equals("+-");
	}
	public String getOpposite()
	{
		if (shortName.equals("<"))
			return ">";
		if (shortName.equals(">"))
			return "<";
		if (shortName.equals("<="))
			return ">=";
		if (shortName.equals(">="))
			return "<=";

		return shortName;
	}

	public String toString()
	{
		return shortName;
	}

}
