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

import java.io.Serializable;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.JoinTable;
import javax.persistence.ManyToMany;
import javax.persistence.ManyToOne;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.hibernate.Criteria;
import org.hibernate.FlushMode;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.annotations.Loggable;
import qspr.util.UserContributedEntity;

import com.eadmet.utils.OCHEMUtils;

/**
 * 
 * A structural alert ("toxicophore", "toxicological alert")
 * In general, its a molecule pattern (SMARTS) related to a property (endpoint) and a publication
 * 
 * @author midnighter
 *
 */

@Entity
@Loggable
@XmlRootElement(name = "substructure-alert")
public class SubstructureAlert implements Serializable, UserContributedEntity
{
	private static final long serialVersionUID = 1L;

	@Id
	@GeneratedValue
	@Column(name = "sa_id")
	@XmlAttribute
	public Long id;

	@Column
	@XmlElement
	public String smart;

	@ManyToOne
	@JoinColumn(name = "property_id")
	public Property property;

	@ManyToOne
	@JoinColumn(name = "article_id")
	@Loggable
	public Article article;

	@XmlTransient
	@Loggable(exclude = true)
	public byte[] image;

	@Column(name = "art_page_num")
	@XmlElement(name = "art-page-num")
	@Loggable(name = "page number")
	public Integer artPageNum;

	@Column(name = "art_table_num")
	@XmlElement(name = "art-table-num")
	@Loggable(name = "table number")
	public String artTableNum;

	@Column(name = "art_line_num")
	@XmlElement(name = "art-line-num")
	@Loggable(name = "line number")
	public Integer artLineNum;

	@Column(name = "art_mol_id")
	@XmlElement(name = "art-mol-id")
	@Loggable(name = "molecule identifier")
	public String artMolId;

	@Column
	@XmlTransient
	public Integer number;

	@Column
	@XmlAttribute
	public boolean folder;

	@Column
	@XmlElement
	public String comment;

	@Column
	@XmlElement
	public String description;

	@Column(name = "smarts_description")
	@XmlElement
	public String smartsDescription;

	@Column
	@XmlElement
	public String name;

	@ManyToOne
	@JoinColumn(name = "modifier_id")
	@XmlTransient
	@Loggable(name = "modifier")
	public User owner; // This is actually _modifier_, not an owner. The old name stayed

	@ManyToOne
	@JoinColumn(name = "introducer_id")
	@XmlTransient
	public User introducer;

	@ManyToOne
	@JoinColumn(name = "con_set_id")
	public ConditionSet conditions;

	@Column(name = "sa_md5", unique = true)
	@XmlAttribute
	@Loggable(exclude = true)
	public String md5;

	@Column(name = "time_modified")
	@XmlTransient
	@Loggable(exclude = true)
	public Timestamp time;

	@Column(name = "time_created")
	@XmlTransient
	@Loggable(exclude = true)
	public Timestamp timeCreated;

	@Column
	@XmlElement
	public Boolean approved;

	@Column
	public int rights;

	@ManyToMany
	(
			targetEntity = SubstructureAlert.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE},
			fetch = FetchType.EAGER
			)
	@JoinTable
	(
			name="SubstructureAlertParent",
			joinColumns={@JoinColumn(name="child_id")},
			inverseJoinColumns={@JoinColumn(name="parent_id")}
			)
	@org.hibernate.annotations.IndexColumn(name="number")
	@XmlTransient
	public List<SubstructureAlert> parents = new ArrayList<SubstructureAlert>();

	@Transient
	public Integer countOccurences;

	/**
	 * How is this alert printed in URLs?
	 */
	@XmlElement
	@Transient
	public String urlFriendlyName;

	@XmlElement(name = "article-abbr")
	protected String getArticleAbbr()
	{
		if (article == null)
			return null;
		Calendar calendar = Calendar.getInstance();
		calendar.setTime(article.publicationDate);
		return "" + (calendar.get(Calendar.YEAR) + 1900) + " " + article.authors.get(0).lastName;
	}

	@XmlElementWrapper(name = "parents")
	@XmlElement(name = "parent")
	protected List<SubstructureAlert> getParents()
	{
		if (parents.isEmpty())
			return null;
		List<SubstructureAlert> uiParents = new ArrayList<SubstructureAlert>();
		for (SubstructureAlert alert : parents)
		{
			SubstructureAlert uiAlert = new SubstructureAlert();
			uiAlert.name = alert.name;
			uiAlert.id = alert.id;
			uiParents.add(uiAlert);
		}

		return uiParents;
	}

//	public boolean isValidSMART()
//	{
//		try
//		{
//			MolImporter.importMol("[$([*R2]([*R])([*R])([*R]))].[$([*R2]([*R])([*R])([*R]))]", QSPRConstants.CXSMARTS);
//			return true;
//		}
//		catch (MolFormatException e)
//		{
//			return false;
//		}
//	}

	/**
	 * Update the internal hash of this SMART used to identify and avoid duplicates
	 */
	public void updateHash()
	{
		StringBuffer hash = new StringBuffer();
		hash.append(smart);
		hash.append(article.id);
		hash.append(property.id);

		if (folder)
			hash.append(name);

		if (rights == Globals.RIGHTS_NONE)
			hash.append(Globals.userSession().user.id);

		md5 = OCHEMUtils.getMD5(hash.toString());
	}

	@XmlElement(name = "hasImage")
	protected boolean getHasImage()
	{
		return image != null;
	}

	public transient SubstructureAlert dublicate;

	/**
	 * Is there a duplicate of this alert?
	 */
	public boolean hasConflicts()
	{
		if (md5 == null)
			return false;

		Globals.session().setFlushMode(FlushMode.MANUAL); //Added to avoid session flush before query to database. Session flush may result in dublicate constraint violation.
		Criteria criteria = Globals.session().createCriteria(SubstructureAlert.class).add(Restrictions.eq("md5", md5));
		if (id != null)
			criteria.add(Restrictions.ne("id", this.id));

		dublicate = null;
		@SuppressWarnings("unchecked")
		List<SubstructureAlert> results = criteria.list();
		if (results.size() > 0)
			dublicate = results.get(0);

		Globals.session().setFlushMode(FlushMode.AUTO);

		return (dublicate != null);
	}

	@XmlAttribute(name = "selected")
	protected boolean isSelected()
	{
		return id != null && Globals.userSession() != null && Globals.userSession().selectedAlerts.contains(id);
	}

	@XmlElement
	protected String getTime()
	{
		return (time != null)?ThreadScope.get().fullDateFormat.format(time):"";
	}

	@XmlElement(name = "time-created")
	protected String getTimeCreated()
	{
		return (timeCreated != null)?ThreadScope.get().fullDateFormat.format(timeCreated):"";
	}

	@XmlElement(name = "owner")
	protected String getOwnerStr()
	{
		if (owner != null)
			return owner.login;
		return null;
	}

	@XmlElement(name = "introducer")
	protected String getIntroducerStr()
	{
		if (introducer != null)
			return introducer.login;
		return null;
	}

	public static SubstructureAlert getByID(Long id)
	{
		return (SubstructureAlert) Globals.session().get(SubstructureAlert.class, id);
	}

	@Override
	public Integer getRights() 
	{
		return rights;
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

	@XmlTransient
	public String getFullSMARTS()
	{
		return AlertSubstitutionVariable.substituteVariables(smart);
	}

	/**
	 * Get a simplified reresentation of an alert ready for marshalling
	 * @return
	 */
	@XmlTransient
	public SubstructureAlert getSimpleAlert()
	{
		SubstructureAlert uiAlert = new SubstructureAlert();
		uiAlert.id = id;
		uiAlert.name = name;
		uiAlert.folder = folder;

		return uiAlert;

	}
}
