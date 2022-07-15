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
import java.util.Calendar;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

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
import javax.persistence.OneToMany;
import javax.persistence.OneToOne;
import javax.persistence.OrderBy;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.FlushMode;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.ProjectionList;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.OCHEMConfiguration;
import qspr.SessionVariable;
import qspr.ThreadScope;
import qspr.annotations.Loggable;
import qspr.dao.Repository;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.applier.ModelApplierTaskProcessor;
import qspr.util.AccessChecker;
import qspr.util.DynaWrap;
import qspr.util.HashedEntity;
import qspr.util.UserContributedEntity;
import qspr.util.ValueDistribution;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;

@Entity
@XmlRootElement(name = "property")
@Loggable(greedy = true)
@SuppressWarnings("unchecked")
public class Property implements UserContributedEntity, HashedEntity<Property>
{

	@Id
	@Column(name="property_id")
	@GeneratedValue
	@XmlAttribute
	public Long id;

	@XmlTransient
	@Transient
	public Property persistentProperty;

	private String name;

	private String shortName;

	@Column
	@XmlTransient
	public String description;

	@ManyToOne
	@JoinColumn(name = "unit_id")
	public Unit defaultUnit;

	@OneToOne
	@JoinColumn(name="category_id")
	@XmlTransient
	public UnitCategory unitCategory;

	@Column(name = "is_condition")
	public boolean isCondition;

	@ManyToOne
	@JoinColumn(name = "modifier_id")
	@XmlTransient
	@Loggable(name = "modifier")
	public User owner;

	@ManyToOne
	@JoinColumn(name = "introducer_id")
	@XmlTransient
	public User introducer;

	@ManyToOne//(fetch = FetchType.LAZY)
	@JoinColumn(name = "moderator_id")
	@XmlTransient
	public User moderator;

	@Column
	@XmlAttribute
	public Integer rights = Globals.RIGHTS_NONE;

	@Column(name = "is_directory")
	public boolean isDirectory;

	@ManyToOne
	@JoinColumn(name = "parent_id")
	public Property parent;

	public Boolean approved = false;

	@Column(name = "property_md5")
	@XmlElement
	public String md5;

	@Column(name = "doc_term")
	@XmlElement
	public String documentationTerm;

	@Column
	@XmlTransient
	public String hash;

	@Transient
	@XmlTransient
	public Property duplicate;

	@ManyToMany
	(
			targetEntity = Property.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE}
			//fetch = FetchType.EAGER			
			)
	@JoinTable
	(
			name="ObligatoryCondition",
			joinColumns={@JoinColumn(name="property_id")},
			inverseJoinColumns={@JoinColumn(name="condition_id")}
			)
	@XmlTransient
	public Set<Property> obligatoryConditions = new HashSet<Property>();

	@OneToMany (fetch = FetchType.LAZY, mappedBy = "property")
	@XmlTransient
	public List<ExperimentalProperty> experimentalProperties;

	@Transient
	@XmlElement
	public ValueDistribution distribution; // XML only

	@Column(name = "bonus_points_weight")
	public double bonusPointsWeight;

	@XmlElement(name = "unitCategory")
	public UnitCategory getUnitCategory()
	{
		if (Globals.getMarshallingOption(MarshallingOption.PROPERTY_UNITCATEGORY))
			return unitCategory;
		if (!ThreadScope.get().controller.equals("batchupload") && !ThreadScope.get().controller.startsWith("model"))
			return unitCategory;

		return null;
	}

	@XmlElementWrapper(name = "predicates")
	@XmlElement(name = "predicate")
	public List<Predicate> getPredicates()
	{
		if (Globals.getMarshallingOption(MarshallingOption.PROPERTY_PREDICATES))
			return Globals.session().createCriteria(Predicate.class).list();
		else
			return null;
	}

	@XmlElement
	public String getDescription()
	{
		if (description == null)
			return null;
		else
			return description.replaceAll("\\\"", "");
	}

	@ManyToMany
	(
			targetEntity = Tag.class,
			cascade={CascadeType.PERSIST, CascadeType.MERGE},
			fetch = FetchType.LAZY
			)
	@JoinTable
	(
			name="PropertyTag",
			joinColumns={@JoinColumn(name="property_id")},
			inverseJoinColumns={@JoinColumn(name="tag_id")}
			)
	@XmlTransient
	public Set<Tag> tags = new HashSet<Tag>();

	@OneToMany(mappedBy = "property", cascade = CascadeType.ALL, fetch = FetchType.LAZY, orphanRemoval = true)
	@XmlTransient
	@OrderBy("name")
	public List<PropertyOption> options = new ArrayList<PropertyOption>();

	@Column
	public String aliases;


	//For PropertiesController/edit.do - list of used conditions
	@XmlAttribute(name="timesUsed")
	@Transient
	public Long count; // For conditions - a number ot times the've been used // For property - times used in a given unit

	@XmlElement
	@Transient
	public Long countExcluded;

	@XmlElement
	@Transient
	public Long countUniqueCompounds;

	@XmlElementWrapper(name="conditionsUsed")
	@XmlElement(name="condition")
	@Transient
	public List<Property> conditionsUsed; //For property - a list of all conditions ever used for this property
	//End

	/**
	 * one of three types 
	 */

	public static final int TYPE_NUMERIC = 0;
	public static final int TYPE_QUALITATIVE = 1;
	public static final int TYPE_TEXTUAL = 2;

	@Column
	@XmlAttribute
	public int type;

	@XmlElement
	@Transient
	public Unit selectedUnit; // for UI only

	@Transient
	@XmlElement
	public int[] transformationsCount;

	@XmlElement(name = "option")
	public List<PropertyOption> getOptions()
	{
		if (Globals.getMarshallingOption(MarshallingOption.PROPERTY_OPTIONS_APPLIER)){
			ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);

			PropertyOptionsFilter filter = new PropertyOptionsFilter();

			if(applier.modelTasks.size() == 1)
				filter.basketId = applier.modelTasks.get(0).model.trainingSet.id;
			else
			{
				filter.basketIds = new HashSet<Long>();
				for(ModelApplierTaskProcessor task: applier.modelTasks)
					filter.basketIds.add(task.model.trainingSet.id);
			}

			filter.propertyId = id;
			return getOptions(filter);
		}


		if (Globals.getMarshallingOption(MarshallingOption.PROPERTY_OPTIONS) || Globals.getMarshallingOption(MarshallingOption.PROPERTY_OPTIONS_FULL))
			if (id != null)
			{
				Criteria c = Globals.session().createCriteria(PropertyOption.class).add(Restrictions.eq("property", this));
				c.addOrder(Order.asc("name"));
				if (!Globals.getMarshallingOption(MarshallingOption.PROPERTY_OPTIONS_FULL))
					c.setMaxResults(100);
				return c.list();

			}
			else
				return options;
		else
			return null;
	}

	static Basket savedBasket = null;

	//Temporary... will move to PropertyOptionsService
	static public List<PropertyOption> getOptions(PropertyOptionsFilter filter)
	{
		if (filter.propertyId == null)
			return new ArrayList<PropertyOption>();

		Property property = (Property) Globals.session().get(Property.class, filter.propertyId);

		if (property == null)
			return new ArrayList<PropertyOption>();

		if (filter.basketId == null && filter.basketIds == null)
			return property.options;

		List<Object[]> list = filter.basketId != null ? getOptions(filter.basketId, property) : getOptions(filter.basketIds, property);

		List<PropertyOption> options = new ArrayList<PropertyOption>();
		Map<Double,PropertyOption> map = new TreeMap <Double,PropertyOption>();

		for (Object[] objects : list) 
		{
			PropertyOption po = (PropertyOption) objects[0];
			po.countOfRecords = (Long) objects[1];
			map.put(-po.countOfRecords - 1./po.id, po);
		}

		for(Map.Entry<Double,PropertyOption> entry : map.entrySet()) {
			PropertyOption po = entry.getValue();
			options.add(po);
			logger.info("Option " + po + " is used in " + po.countOfRecords + " records");
		}

		return options;
	}

	static private List<Object[]>  getOptions(Set<Long> basketIds, Property property){

		List<Object[]>  list = null;

		for(Long id:basketIds) {

			List<Object[]>  newlist = getOptions(id, property);

			if(list == null) list = newlist;
			else
				for (Object[] newobjects : newlist) 
				{
					boolean found = false;
					PropertyOption npo = (PropertyOption) newobjects[0];

					for (Object[] objects : list) 
					{
						PropertyOption po = (PropertyOption) objects[0];
						if(npo.id != po.id) continue;
						found = true;
						objects[1] = (Long) objects[1] + (Long) newobjects[1];
					}

					if(found) continue;
					list.add(newobjects);
				}
		}

		return list;		
	}


	static private List<Object[]>  getOptions(Long basketId, Property property){

		//Basket basket = Basket.getBasket(Globals.userSession(), basketId);
		// TODO - is this security problem - we do not check permissions to get options of the private basket.... I do not think so
		Basket basket = Repository.basket.getById(basketId);

		List<Object[]>  list= Globals.session().createCriteria(BasketEntry.class)
				.createAlias("ep", "ep1")
				.createAlias("ep1.conditions", "cs")
				.createAlias("cs.values", "pv")
				.add(Restrictions.eq("pv.property", property))
				.add(Restrictions.isNotNull("pv.option"))
				.add(Restrictions.eq("basket", basket))
				.setProjection(Projections.projectionList()
						.add(Projections.groupProperty("pv.option"))
						.add(Projections.countDistinct("ep1.id")))
				.list();

		return list;
	}


	/**
	 * Get the number of options for a qualitative property
	 * @return
	 */
	@XmlElement(name = "options-count")
	public long getOptionsCount()
	{
		if (id == null)
			return 0;
		return (Long) Globals.session().createCriteria(PropertyOption.class).add(Restrictions.eq("property", this)).setProjection(Projections.count("id")).uniqueResult();
	}

	@XmlAttribute(name="property-record")
	public Long getPropertyCount()
	{
		if (Globals.getMarshallingOption(MarshallingOption.PROPERTY_RECORD_COUNT))
		{
			if(id != null)
			{
				long time = Calendar.getInstance().getTimeInMillis();
				Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class)
						.add(Restrictions.isNull("deleted"));
				ExperimentalProperty.addAccessRestrictions(criteria, Globals.RIGHTS_FREELY_AVAILABLE, null, false, true);

				if (isCondition == false)
				{
					if (isDirectory)
						criteria.createCriteria("property").add(Restrictions.eq("parent", this));
					else
						criteria.add(Restrictions.eq("property", this));
				}
				else
					criteria.createCriteria("conditions").createCriteria("values").add(Restrictions.eq("property", this));

				criteria.setProjection(Projections.countDistinct("id"));
				logger.info(name + ": counted " + criteria.list().get(0) + " records, " + (Calendar.getInstance().getTimeInMillis() - time) + " ms.");
				return (Long)criteria.list().get(0);
			}
		}
		return 0L;
	}

	public PropertyOption getOption(String name)
	{
		return (PropertyOption) Globals.session().createCriteria(PropertyOption.class).add(Restrictions.eq("property", this)).add(Restrictions.like("name", name)).uniqueResult();
	}

	@XmlAttribute
	public Boolean isReferenced()
	{
		//if (isCondition == null)
		//	isCondition = false;
		if (Globals.getMarshallingOption(MarshallingOption.PROPERTY_RECORD_COUNT) && id != null && isCondition)
		{
			Long count = (Long) Globals.session()
					.createQuery("select count(*) from ExperimentalProperty ep join ep.conditions cs join cs.values vs where vs.property=:property")
					.setParameter("property", this)
					.uniqueResult();
			return count > 0;
		}
		return null;
	}

	@XmlElement(name = "used-unit")
	public List<Unit> getUsedUnit()
	{
		if (Globals.getMarshallingOption(MarshallingOption.PROPERTY_UNITS))
		{
			if (id == null)
				return null;

			List<Unit> units = new ArrayList<Unit>();

			ProjectionList projList = Projections.projectionList();
			projList.add(Projections.groupProperty("unit"));
			projList.add(Projections.countDistinct("id"));

			Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class)
					.add(Restrictions.eq("property", this))
					.setProjection(projList);

			ExperimentalProperty.addAccessRestrictions(criteria, Globals.RIGHTS_FREELY_AVAILABLE, Globals.userSession().user, false, true);
			List<Object[]> res = criteria.list();
			for (Object[] objects : res) 
			{
				Unit unit = (Unit) objects[0];
				if (unit != null)
				{
					Unit copy = new Unit();
					copy.setName(unit.getName());
					copy.id = unit.id;
					copy.countInProperty = (Long) objects[1];
					units.add(copy);
				}
			}

			return units;
		}

		return null;
	}

	public PropertyOption getOptionByName(String name)
	{
		for (PropertyOption co : options)
			if (co.name.equalsIgnoreCase(name))
				return co;
		return null;
	}

	public static Property getById(Long id)
	{
		return (Property) Globals.session().get(Property.class, id);
	}

	public static Property getByName(String name)
	{
		Criteria c = Globals.session().createCriteria(Property.class).add(Restrictions.eq("shortName", Property.shortName(name)));
		List<Property> p = c.list();
		if (p.size() > 0)
			return p.get(0);
		return null;
	}

	public Property()
	{
	}

	public Property(String _name)
	{
		setName(_name);
	}

	public String toString()
	{
		return this.name;
	}

	@XmlAttribute@Column(unique = true)
	public Property setName(String _name)
	{
		name = _name.trim();
		shortName = Property.shortName(name);
		updateHash();
		return this;
	}

	public String getName()
	{
		return name;
	}

	public String getShortName()
	{
		return shortName;
	}

	@XmlElement(name = "owner")
	protected String getOwnerStr()
	{
		if (owner != null)
			return owner.login;
		else
			return null;
	}

	@XmlElement(name = "introducer")
	protected String getIntroducerStr()
	{
		if (introducer != null)
			return introducer.login;
		else
			return null;
	}

	@XmlElement(name = "moderator")
	protected User getModeratorXml()
	{
		if (Globals.getMarshallingOption(MarshallingOption.PROPERTY_MODERATOR))
			return moderator;
		return null;

	}

	@Override
	public boolean equals(Object obj)
	{
		if (obj != null && !(obj instanceof Property))
			return false;
		Property property = (Property) obj;
		if (this.id != null)
			return property != null && this.id.equals(property.id);
		else
			return this == obj;
	}

	@Override
	public int hashCode()
	{
		if (this.id != null)
			return this.id.hashCode();
		else
			return super.hashCode();
	}

	public static String shortName(String name)
	{
		String result = name.toLowerCase().replaceAll("[^0-9a-z\\+\\-\\%]*", "");
		return result;
	}

	public boolean isNumeric()
	{
		return type == TYPE_NUMERIC;
	}


	public boolean isQualitative()
	{
		return type == TYPE_QUALITATIVE;
	}

	public boolean isTextual()
	{
		return type == TYPE_TEXTUAL;
	}

	@XmlAttribute
	protected String getQualitive()
	{
		// For XML backward-compatibility
		if (isQualitative())
			return "true";
		return "false";
	}

	@XmlAttribute(name = "children-count")
	protected Long getChildrenCount()
	{
		if (!isDirectory)
			return 0L;
		if (id == null)
			return 0L;
		return (Long) Globals.session().createCriteria(Property.class).add(Restrictions.eq("parent", Globals.session().get(Property.class, this.id))).setProjection(Projections.countDistinct("id")).uniqueResult();
	}

	public void assignModerator(User user)
	{
		// Notify the newly-chosen moderator
		if (user.isExtended()) {
			DynaWrap extended = user.dynaWrapped();
			Mailer.postMailSafely(new Email(extended.getString("email"), "You are the moderator of \"" + name + "\" at ochem.eu", "Dear " + extended.getString("title") + " " + extended.getString("firstName") + " " + extended.getString("lastName") + ",\n\n"
					+"you have been assigned as a moderator of all the data for \""+name+"\" on OCHEM resource.\n"+
					"You can approve or disapprove new data published by other users for this property.\n\nBest regards,\nOCHEM Team"));
		}

		// Notify administrators
		Mailer.notifyAdmins(user.login + " assigned a moderator of a property " + name, "We thought it might be interesting for you that user " + user.login + " has been assigned as a moderator of a property " + name);

		moderator = user;
	}

	@Override
	public Integer getRights() {
		return rights;
	}

	@Override
	public User getIntroducer() {
		return introducer;
	}

	@Override
	public User getOwner() {
		return owner;
	}

	@Override
	public void updateHash() 
	{
		hash = shortName;
		//		if (rights == Globals.RIGHTS_NONE)
		//			if (owner.group != null)
		//				hash += "_g"+owner.group.id;
		//			else
		//				hash += "_u"+owner.id;
		md5 = OCHEMUtils.getMD5(hash);
	}

	@Override
	public boolean hasConflicts() 
	{
		if (md5 == null)
			return false;

		Globals.session().setFlushMode(FlushMode.MANUAL); //Added to avoid session flush before query to database. Session flush may result in dublicate constraint violation.
		Criteria criteria = Globals.session().createCriteria(Property.class).add(Restrictions.eq("md5", md5));
		if (id != null)
			criteria.add(Restrictions.ne("id", this.id));

		duplicate = null;
		List<Property> eps = criteria.list();
		if (eps.size() > 0)
			duplicate = eps.get(0);

		Globals.session().setFlushMode(FlushMode.AUTO);

		return (duplicate != null);
	}

	@Override
	public Property getDuplicate() 
	{
		return duplicate;
	}

	@XmlTransient
	public boolean isPublished()
	{
		return rights == Globals.RIGHTS_FREELY_AVAILABLE && approved != null && approved;
	}

	public void approve()
	{
		AccessChecker.requestModeratorPrivileges();
		approved = true;
	}

	public void unapprove()
	{
		AccessChecker.requestModeratorPrivileges();
		approved = false;
		if (id != null)
		{
			Globals.session().createSQLQuery("update ExperimentalProperty set approved = false where property_id=" + id).executeUpdate();
			Globals.session().saveOrUpdate(this);
		}

	}

	public void publish()
	{
		rights = Globals.RIGHTS_FREELY_AVAILABLE;
	}

	static public long getDummyId(){
		if(DUMMYID == null)
			DUMMYID= Repository.property.getProperty(QSPRConstants.DUMMY, false).id;
		if(DUMMYID == null) throw new CriticalException("The database should contain property " + QSPRConstants.DUMMY);
		return DUMMYID;
	}

	private static transient final Logger logger = LogManager.getLogger(Property.class);

	public static Long DUMMYID = null;

	/**
	 * So far there is only one mappable condition
	 * Can be extended later to all of them, such as Catalizer, Ligand, Reagent, etc.
	 * @return article which is used to map name, SMILES to canonical name
	 */
	
	public Long isMapable() {
		if(name.equalsIgnoreCase(QSPRConstants.SOLVENT_CONDITION))
			return OCHEMConfiguration.solvent;
		return null;
	}
}
