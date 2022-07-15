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

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.ManyToMany;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;

@Entity
@XmlRootElement(name = "author")
public class Author 
{

	@Id
	@Column(name = "author_id")
	@GeneratedValue
	@XmlAttribute
	public Long id;


	@Column(name = "first_name")
	@XmlElement(name = "FirstName")
	public String firstName;

	@Column(name = "last_name")
	@XmlElement(name = "LastName")
	public String lastName;

	@Column(name = "initials")
	@XmlElement(name = "Initials")
	public String initials;

	@XmlAttribute(name = "printed-name")
	public String getPrintedName() 
	{
		if (id == null)
			return "New author*";
		else
		{
			if((initials != null) && (!initials.equals("")))
				return lastName+", "+initials;
			else
				return lastName;
		}

	}

	@ManyToMany(mappedBy = "authors")
	Set<Article> articles = new HashSet<Article>();

	public static Author get(String firstName, String lastName, String initials)
	{
		Criteria criteria = Globals.session().createCriteria(Author.class).add(Restrictions.eq("lastName", lastName));

		if (!firstName.equals(""))
			criteria.add(Restrictions.eq("firstName", firstName));

		if (!initials.equals(""))
			criteria.add(Restrictions.eq("initials", initials));

		@SuppressWarnings("unchecked")
		List<Author> authors = (List<Author>)criteria.list();

		if (authors.size() > 0)
			return authors.get(0);
		else
		{
			Author author = new Author();
			author.firstName = firstName;
			author.lastName = lastName;
			author.initials = initials;
			Globals.session().saveOrUpdate(author);
			return author;
		}
	}

	public static Author get(String title)
	{
		String[] parts;
		if(title.contains(","))
			parts = title.split(",");
		else
			parts = title.split(" ");
		return getAuthorName(parts);
	}

	public static Author getFromRIS(String title)
	{
		String[] parts = title.split(",");
		return getAuthorName(parts);		
	}

	public static List<Author> getFromISI(String value)
	{
		String[] authors = value.split("#@#");
		List<Author> res = new ArrayList<Author>();
		for (int i=0; i<authors.length; i++)
		{
			if (!authors[i].trim().equals(""))
			{
				String[] pieces = authors[i].split(",");
				res.add(getAuthorName(pieces));
			}
		}
		return res;
	}

	public static List<Author> getFromEN(String value)
	{
		List<Author> res = new ArrayList<Author>();
		if(value.contains(" "))
			res.add(getAuthorName(value.split(" ")));
		else if(value.contains(","))
			res.add(getAuthorName(value.split(",")));
		return res;
	}

	public static List<Author> getFromBT(String value)
	{
		String[] authors = value.split(" and ");
		List<Author> res = new ArrayList<Author>();
		for (int i=0; i<authors.length; i++)
		{
			if (!authors[i].trim().equals(""))
			{
				String[] pieces = authors[i].split(",");
				String _lastName = pieces[0];
				if (!_lastName.equals("others"))
				{
					res.add(getAuthorName(pieces));
				}
			}
		}
		return res;
	}

	private static Author getAuthorName(String[] parts) {
		String _lastName = parts[0].trim();
		String _firstName = "";
		String _initials = "";
		if (parts.length > 2)
		{
			for(int i=1; i < parts.length; i++)
			{
				if(parts[i].length() > 2)
					_firstName = parts[i];
				else
					_initials = parts[i];
			}
		}
		else if (parts.length == 2)
		{
			String string =  parts[parts.length - 1].trim();
			if(string.length() > 2)
			{
				_firstName = string;
				_initials = string.substring(0, 1)+".";
			}
			else	 
				_initials = string;
		}
		return Author.get(_firstName, _lastName, _initials);
	}

	public String toString()
	{
		return this.lastName;
	}

	static Author getByName(String firstName, String lastName, String initials, boolean createIfMissing)
	{
		@SuppressWarnings("unchecked")
		List<Author> l = Globals.session().createCriteria(Author.class)
		.add(Restrictions.like("firstName",firstName))
		.add(Restrictions.like("lastName",lastName))
		.add(Restrictions.like("initials",initials))
		.list();
		if (l.size() > 0)
			return l.get(0);

		if (!createIfMissing)
			return null;

		Author au = new Author();
		au.firstName = firstName;
		au.lastName = lastName;
		au.initials = initials;
		return au;
	}

}
