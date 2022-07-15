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

import java.util.Date;
import java.util.List;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.annotations.Loggable;

import com.eadmet.exceptions.UserFriendlyException;

@Entity
@XmlRootElement(name = "journal")
//@Loggable
public class Journal 
{
	private static transient final Logger logger = LogManager.getLogger(Journal.class);

	@Id
	@Column(name = "journal_id")
	@GeneratedValue
	@XmlAttribute
	public Long id;

	@XmlTransient
	private String title;

	@XmlElement
	public String publisher;

	@XmlTransient
	public Date publish_date;

	@XmlTransient
	private String abbreviation;

	@XmlTransient
	private String issn;

	@ManyToOne
	@JoinColumn(name = "modifier_id")
	@XmlTransient
	@Loggable(name = "modifier")
	public User owner;

	@ManyToOne
	@JoinColumn(name = "introducer_id")
	@XmlTransient
	public User introducer;

	@XmlElement
	private String link;

	@XmlAttribute(name="articlesCount")
	public Long getReferencesCount()
	{
		if (ThreadScope.get().controller.equals("journal") && this.id != null)
		{
			Criteria criteria = Globals.session().createCriteria(Article.class)
					.add(Restrictions.eq("journal", this));
			criteria.setProjection(Projections.rowCount());
			return (Long)criteria.list().get(0);
		}
		return 0L;
	}

	@XmlElement(name = "owner")
	public String getOwner()
	{
		if (owner != null)
			return owner.login;
		else
			return null;
	}

	@XmlElement(name = "introducer")
	public String getIntroducer()
	{
		if (introducer != null)
			return introducer.login;
		else
			return null;
	}

	private static String preprocessISSN(String issn)
	{
		if (issn != null && !issn.equals(""))
			return issn.replace("-", "").replace(" ", "").trim();
		else
			return null;
	}

	private static String preprocessTitle(String title)
	{
		if (title != null && !title.equals(""))
			return title.trim();
		else
			return null;
	}

	@XmlElement(name = "abbreviation")
	public String getAbbreviation()
	{
		return abbreviation;
	}

	public void setAbbreviation(String abbr)
	{
		if (abbr != null && !abbr.equals(""))
			abbreviation = abbr.trim();
		else
			abbreviation = null;
	}


	@XmlElement(name = "title")
	public String getTitle()
	{
		return title;
	}

	public void setTitle(String title)
	{
		this.title = preprocessTitle(title);
	}


	@XmlElement(name = "issn")
	public String getISSN()
	{
		if (issn != null)
			return issn.substring(0, 4) + "-" + issn.substring(4);
		return null;
	}

	public void setISSN(String issn)
	{
		this.issn = preprocessISSN(issn);
	}


	public void setLink(String link)
	{
		if(!link.equals("") && !link.startsWith("http://"))
			this.link = "http://" + link;
		else
			this.link = link;
	}

	@SuppressWarnings("unchecked")
	public static Journal getByISSN(String issn) throws Exception
	{
		issn = preprocessISSN(issn);
		List<Journal> journals = 
				Globals.session().createCriteria(Journal.class)
				.add(Restrictions.eq("issn", issn))
				.list();

		if (journals.size() > 0)
			return journals.get(0);
		else
			return null;
	}

	@SuppressWarnings("unchecked")
	public static Journal getByTitle(String title) 
	{
		List<Journal> journals = 
				Globals.session().createCriteria(Journal.class)
				.add(Restrictions.eq("title", preprocessTitle(title)))
				.list();

		if (journals.size() > 0)
			return journals.get(0);

		journals = 
				Globals.session().createCriteria(Journal.class)
				.add(Restrictions.eq("abbreviation", preprocessTitle(title)))
				.list();

		if (journals.size() > 0)
			return journals.get(0);

		return null;

	}

	public void doCheckRights() throws Exception
	{	
		if ( ! (this.owner == null || this.owner.equals(Globals.userSession().user)))
		{
			if ( ! (Globals.userSession().user != null && Globals.userSession().user.rank > this.owner.rank))
			{
				if ( ! (Globals.userSession().user != null && Globals.userSession().user.rank.equals(this.owner.rank)))
				{
					logger.info("Not permitted, edit/delete has been ignored");
					throw new UserFriendlyException("Not permitted, edit/delete has been ignored");
				}
			}
		}
	}

	@XmlElement(name = "publish_date")
	XmlDate getPublicationDate()
	{
		return new XmlDate(publish_date);
	}

	public static final String defaultDateFormat = "yyyy-MM-dd HH:mm:ss";

	@SuppressWarnings("unchecked")
	public static Journal getByTitle(String title, boolean createIfMissing)
	{
		List<Journal> journals = Globals.session().createCriteria(Journal.class)
				.add(Restrictions.or(Restrictions.eq("title", title), Restrictions.eq("abbreviation", title)))
				.list();

		if(journals.size() > 0)
			return journals.get(0);

		if (!createIfMissing)
			return null;

		Journal journal = new Journal();
		journal.title = title;
		return journal;			

	}

	public void update(Journal fetchedJournal) {
		if(fetchedJournal == null) return;

		if(fetchedJournal.issn == null || fetchedJournal.issn.length() == 0)
			throw new UserFriendlyException("The journal was not found in NCBI and thus no update will be performed");

		if(issn == null || issn.length() == 0) issn = fetchedJournal.issn;

		if(!issn.equals(fetchedJournal.issn))
			throw new UserFriendlyException("The journal found in NCBI has a diffenent issn (" + 
					fetchedJournal.issn + ") than the query: " + issn);

		title = fetchedJournal.title;
		abbreviation = fetchedJournal.abbreviation;
		publisher = fetchedJournal.publisher;
		publish_date = fetchedJournal.publish_date;
	}

}
