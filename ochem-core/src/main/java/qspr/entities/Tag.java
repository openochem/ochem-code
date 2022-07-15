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

import java.math.BigInteger;
import java.util.Calendar;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToMany;
import javax.persistence.ManyToOne;
import javax.persistence.Transient;
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
import qspr.util.UserContributedEntity;

import com.eadmet.exceptions.UserFriendlyException;

/**
 * A tag for properties or molecules
 * @author midnighter
 *
 */
@Entity
@XmlRootElement(name = "tag")
public class Tag implements UserContributedEntity
{
	private static transient final Logger logger = LogManager.getLogger(Tag.class);

	@Id
	@Column(name = "tag_id")
	@GeneratedValue
	@XmlAttribute
	public Long id;

	@Column
	@XmlAttribute
	public String name;

	@Column
	@XmlAttribute
	public String type;

	@Column
	@XmlElement
	public String description;

	@ManyToOne
	@JoinColumn(name = "modifier_id")
	@XmlTransient
	@Loggable(name = "modifier")
	public User owner;

	@XmlAttribute(name="basketCount")
	@Transient
	public Long count;

	@ManyToOne
	@JoinColumn(name = "introducer_id")
	@XmlTransient
	@Loggable
	public User introducer;

	@Column(name = "doc_term")
	@XmlElement
	public String documentationTerm;

	@Column(name = "is_public")
	@XmlAttribute
	public boolean isPublic;

	@Column(name = "show_in_browser")
	@XmlAttribute
	public boolean showInBrowser;


	@ManyToMany(
			mappedBy = "tags",
			targetEntity = Property.class,
			fetch = FetchType.LAZY
			)
	@XmlTransient
	public Set<Property> properties = new HashSet<Property>();

	@ManyToMany(
			mappedBy = "tags",
			targetEntity = Mapping1.class,
			fetch = FetchType.LAZY
			)
	@XmlTransient
	public Set<Mapping1> mapping = new HashSet<Mapping1>();

	@ManyToMany(
			mappedBy = "tags",
			targetEntity = Article.class,
			fetch = FetchType.LAZY
			)
	@XmlTransient
	public Set<Article> articles = new HashSet<Article>();

	public static Tag getByName(String name)
	{
		Tag tag = (Tag) Globals.session().createCriteria(Tag.class)
				.add(Restrictions.like("name", name)).uniqueResult();

		return tag;
	}

	public static Tag getByID(Long id)
	{
		Tag tag = (Tag) Globals.session().get(Tag.class, id);

		return tag;
	}

	public void doCheckRights() throws Exception
	{
		if (!(this.owner == null || this.owner.equals(Globals.userSession().user)))
			if (!(Globals.userSession().user != null && Globals.userSession().user.rank > this.owner.rank))
			{
				if (!(Globals.userSession().user != null && Globals.userSession().user.rank.equals(this.owner.rank)))
				{
					logger.info("Not permitted, edit has been ignored");
					throw new UserFriendlyException("Not permitted, edit has been ignored");
				}
			}
	}

	public String toString()
	{
		return name;
	}

	public boolean equals(Object obj)
	{

		Tag tag = (Tag) obj;
		return (tag != null) && this.id.equals(tag.id);
	}

	@XmlElement(name = "owner")
	public String getOwnerLogin()
	{
		if (owner != null)
			return owner.login;
		else
			return null;
	}

	@XmlElement(name = "introducer")
	public String getIntroducerLogin()
	{
		if (introducer != null)
			return introducer.login;
		else
			return null;
	}

	public boolean canEdit()
	{
		if (Globals.userSession().user == null)
			return false;

		if (isPublic)
			return true;

		if (Globals.userSession().user.equals(introducer))
			return true;

		if (introducer.group != null)
			return introducer.group.equals(Globals.userSession().user.group);

		return false;
	}

	@XmlTransient
	public int getMoleculesCount()
	{
		if (id == null)
			return 0;
		return ((BigInteger) Globals.session().createSQLQuery("select count(*) from MoleculeTag where tag_id=:id").setLong("id", id).uniqueResult()).intValue();
	}

	@XmlAttribute(name="tag-molecule")
	public Long getRecordSize(){
		if (ThreadScope.get().controller.equals("tags"))
		{
			long time = Calendar.getInstance().getTimeInMillis();
			try
			{
				if((this.id != null) && (this.id > 0)){
					Criteria criteria = Globals.session().createCriteria(Mapping1.class);

					if(this.type.equals("property"))
					{
						@SuppressWarnings("unchecked")
						List<Property> pls = Globals.session().createCriteria(Property.class)
						.createAlias("tags", "t")
						.add(Restrictions.eq("t.id", this.id)).list();

						if (pls.size() == 0)
							return 0L;

						criteria
						.createAlias("molecules", "m")
						.createAlias("m.experimentalProperties", "ep")
						.add(Restrictions.in("ep.property", pls.toArray()));
					}else{
						criteria.createAlias("tags", "t")
						.add(Restrictions.eq("t.id", this.id));
					}
					criteria.setProjection(Projections.countDistinct("id"));
					return (Long)criteria.list().get(0);
				}
			}
			finally
			{
				logger.info("Counting for tag " + name + ": " + (Calendar.getInstance().getTimeInMillis() - time) + "ms.");
			}
		}
		return 0L;
	}

	@Override
	public User getIntroducer() 
	{
		return introducer;
	}

	@Override
	public User getOwner() 
	{
		return owner;
	}

	@Override
	public Integer getRights() 
	{
		if (isPublic)
			return Globals.RIGHTS_FREELY_AVAILABLE;
		else
			return Globals.RIGHTS_NONE;
	}

}
