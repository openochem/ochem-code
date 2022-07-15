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

import java.io.IOException;
import java.io.Serializable;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeoutException;

import javax.persistence.CascadeType;
import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.FetchType;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.persistence.JoinColumn;
import javax.persistence.ManyToOne;
import javax.persistence.OneToMany;
import javax.persistence.Transient;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.Hibernate;
import org.hibernate.SQLQuery;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.type.LongType;

import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.mailer.Mailer;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.business.BasketFilter;
import qspr.business.BasketPeer;
import qspr.business.ModelsFilter;
import qspr.business.Privileges;
import qspr.dao.Repository;
import qspr.util.MoleculePeer;
import qspr.util.UploadContext;
import qspr.util.UserContributedEntity;

@Entity
@XmlRootElement(name = "basket")
@SuppressWarnings("unchecked")
public class Basket implements Serializable, UserContributedEntity
{
	private static transient final Logger logger = LogManager.getLogger(Basket.class);
	private static final long serialVersionUID = 1L;

	@Id
	@GeneratedValue
	@Column(name = "basket_id")
	@XmlAttribute
	public Long id;

	@ManyToOne
	@JoinColumn(name = "user_id")
	@XmlTransient
	public User user;

	@ManyToOne
	@JoinColumn(name = "session_id")
	@XmlTransient
	public Session session;

	@Column
	@XmlAttribute
	public String name;

	@Column
	@XmlElement
	public String description;

	@Column
	@XmlElement
	public String excludedrecords;

	@Column(name = "is_molecule_set")
	@XmlTransient
	public boolean isMoleculeSet;

	@OneToMany(mappedBy = "basket", cascade={CascadeType.ALL}, fetch = FetchType.LAZY)
	@XmlTransient
	public List<BasketEntry> entries = new ArrayList<BasketEntry>();

	@Column(name = "basket_type", nullable = false)
	@XmlTransient
	public Long basketType = 0L; //0 = normal, 1 = system

	@Column
	@XmlAttribute
	public Integer rights = Globals.RIGHTS_NONE; //
	//	@XmlAttribute(name="size")
	//	public Integer getEntrieSize()
	//	{
	//		return entries.size();
	//	}

	//count of property, article and condition to show statistics 
	@XmlElementWrapper(name="propertyUsed")
	@XmlElement(name="property")
	@Transient
	public List<Property> propertyUsed; //For property

	@XmlElementWrapper(name="articleUsed")
	@XmlElement(name="article")
	@Transient
	public List<Article> articleUsed; //For article

	@XmlElementWrapper(name="tagUsed")
	@XmlElement(name="tag")
	@Transient
	public List<Tag> tagUsed; //For tags; currently disabled in BasketService

	@XmlElement
	@Transient
	public Long totalUniqueCompounds_SC; //for basket, total number of unique compounds with stereochemistry 

	@XmlElement
	@Transient
	public Long totalUniqueCompounds; // for basket, total number of unique compounds without stereochemistry

	@XmlTransient
	@Column(name = "cached_count")
	public Long cachedCount;

	@XmlTransient
	@Column(name = "last_modified")
	public Timestamp lastModified;

	@Transient
	public static Basket getBasket(Session session, String name)
	{
		return getBasket(session, name, true);
	}

	public static Basket getById(Long id)
	{
		return (Basket) Globals.session().get(Basket.class, id);
	}

	@Transient
	public static Basket getBasket(Session session, String name, boolean group)
	{
		Criteria criteria = Globals.session().createCriteria(Basket.class);
		criteria.add(Restrictions.eq("name", name));
		if (session.user != null)
			if (session.user.group == null || !group)
				criteria.add(Restrictions.eq("user", session.user));
			else
				// Allow to see baskets of my group
				criteria.createCriteria("user").add(Restrictions.eq("group", session.user.group));
		else
			criteria.add(Restrictions.eq("session", session));
		List<Basket> baskets = criteria.list();
		Basket basket;

		if (baskets.size() == 0)
		{
			basket = new Basket();
			basket.name = OCHEMUtils.getFilteredBasketName(name);
			basket.session = session;
			basket.user = session.user;
			Globals.session().save(basket);
		}
		else
			basket = baskets.get(0);

		return basket;
	}

	@XmlTransient
	public Privileges getPrivileges()
	{
		Privileges privileges = new Privileges("basket");

		privileges.canEdit = OCHEMConfiguration.autoLoginUser != null ? true:
			Globals.userSession().equals(session) || (Globals.userSession().user != null && Globals.userSession().user.equals(session.user));
		privileges.canView = privileges.canEdit;

		if (!privileges.canView)
		{
			// Consider if this basket is "public"
			BasketFilter filter = new BasketFilter();
			filter.id = id;
			filter.showPublicBaskets = filter.showGroupBaskets = true;
			privileges.canView = BasketPeer.getListCriteria(filter).list().size() > 0;
		}

		return privileges;
	}

	@Transient
	public static Basket getBasket(Session session, Long id)
	{
		Basket basket = Repository.basket.getById(id);
		if(basket != null && basket.user != null && basket.user.isPublisher()) return basket;

		Criteria criteria = Globals.session().createCriteria(Basket.class, "basket");
		Disjunction disjunction = Restrictions.disjunction();
		criteria.add(Restrictions.eq("id", id));
		if (session.user != null)
			if (session.user.group == null)
				disjunction.add(Restrictions.eq("user", session.user));
			else
			{
				criteria.createAlias("user", "u");
				disjunction.add(Restrictions.eq("u.group", session.user.group));
			}
		else
			disjunction.add(Restrictions.eq("session", session));
		/*
		// Reveal also baskets, associated with published models / Midnighter
		DetachedCriteria modelCriteria = DetachedCriteria.forClass(Model.class, "model");
		modelCriteria.setProjection(Projections.id());
		modelCriteria.add(Restrictions.or(
				Restrictions.eqProperty("trainingSet.id", "basket.id"),
				Restrictions.eqProperty("validationSet.id", "basket.id")));

		Disjunction modelIsAvailable = Restrictions.disjunction();
		modelIsAvailable.add(Restrictions.eq("model.published", new Boolean(true)));
		if (!session.visitedModelsIds.isEmpty())
			modelIsAvailable.add(Restrictions.in("model.publicId", session.visitedModelsIds));
		modelCriteria.add(modelIsAvailable);

		disjunction.add(Subqueries.exists(modelCriteria));
		 */		

		if (session.user == null || !session.user.isSuperUser())
			criteria.add(disjunction);
		List<Basket> baskets = criteria.list();

		return baskets.size() == 0 ? null : baskets.get(0);
	}

	public static Basket getBasket(String name)
	{
		return getBasket(name, Globals.userSession(), true);
	}

	public static Basket getBasket(String name, Session session)
	{
		return Basket.getBasket(name, session, true);
	}

	public static Basket getBasket(String name, Session session, boolean createIfNotExists)
	{
		Criteria c = Globals.session()
				.createCriteria(Basket.class)
				.add(Restrictions.eq("name", name));

		if (session.user != null)
			c.add(Restrictions.eq("user", session.user));
		else
			c.add(Restrictions.eq("session", session));

		List<Basket> baskets = c.list(); 
		Basket basket;
		if (baskets.size() == 0)
		{
			if (!createIfNotExists)
				return null;
			logger.info("Creating a new basket with name " + name);
			basket = new Basket();
			basket.name = OCHEMUtils.getFilteredBasketName(name);
			basket.session = session;
			basket.user = basket.session.user;
			Globals.session().saveOrUpdate(basket);
		}
		else
		{
			basket = baskets.get(0);
			basket.session = session;
			Globals.session().saveOrUpdate(basket);
		}

		return basket;
	}

	public static Basket getBasket(String name, UploadContext context) throws Exception
	{
		Basket basket = null;
		try
		{
			basket = context.basketCache.get(name);

			if (basket != null)
				return basket;

			basket = getBasket(name);
			return basket;

		} catch (Exception e)
		{
			context.basketCache.put(name,e);
			throw e;
		} finally
		{
			if (basket != null)
				context.basketCache.put(name, basket);
		}
	}

	public void countRecordsByProperties(boolean excludedRecordsOnly)
	{
		// Enumerate all the properties in the basket
		Criteria propertyCriteria = Globals.session().createCriteria(BasketEntry.class);
		propertyCriteria.add(Restrictions.eq("basket", this));
		propertyCriteria.createAlias("ep", "e");
		propertyCriteria.createAlias("e.molecule", "m");
		propertyCriteria.createAlias("m.mapping2", "mp2");
		propertyCriteria.setProjection(Projections.projectionList()
				.add(Projections.groupProperty("e.property"))
				.add(Projections.countDistinct("id"))
				.add(Projections.countDistinct("mp2.id")));

		if (excludedRecordsOnly)
			propertyCriteria.add(Restrictions.eq("exclude", Boolean.TRUE));

		List<Object[]> propertyList = propertyCriteria.list();
		if (propertyList.isEmpty())
			return;

		if (this.propertyUsed == null)
			this.propertyUsed = new ArrayList<Property>();

		// Count records by properties
		for (Object[] objects : propertyList)
		{
			Property property = (Property) objects[0];
			if (!excludedRecordsOnly)
			{
				property.count = (Long)objects[1]; 
				property.countUniqueCompounds = (Long) objects[2];
			}
			else
				property.countExcluded = (Long) objects[1];
			if (!this.propertyUsed.contains(property))
				this.propertyUsed.add(property);

			if (property.isQualitative())
			{
				// Count individual options of a property
				propertyCriteria = Globals.session().createCriteria(BasketEntry.class);

				propertyCriteria.createAlias("ep", "e");
				propertyCriteria.createAlias("e.molecule", "m");
				propertyCriteria.createAlias("m.mapping2", "mp2");

				propertyCriteria.add(Restrictions.eq("basket", this));
				propertyCriteria.add(Restrictions.eq("e.property", property));

				propertyCriteria.setProjection(Projections.projectionList()
						.add(Projections.groupProperty("e.option"))
						.add(Projections.countDistinct("id"))
						.add(Projections.countDistinct("mp2.id")));

				if (excludedRecordsOnly)
					propertyCriteria.add(Restrictions.eq("exclude", Boolean.TRUE));

				List<Object[]> optionsList = propertyCriteria.list();
				for (Object[] row : optionsList)
				{
					logger.info(property.getName() + "  " + row[0] + " records: " + row[1] + " unique compounds: " + row[2]);
					PropertyOption option = (PropertyOption) row[0];

					if(option == null) {
						Mailer.notifyDevelopers("Error in property" + property.getName(), property.getName()+ " " + property.id+ " has null option");
						continue;
					}

					if (!excludedRecordsOnly)
					{
						option.countOfRecords = (Long) row[1];
						option.countUniqueCompounds = (Long) row[2];
					}
					else
						option.countExcluded = (Long) row[1];
				}

			}
		}
	}

	public transient List<Property> modelProperties;

	public boolean isEditableBy(Session requestingSession)
	{
		if(OCHEMConfiguration.autoLoginUser != null) return true;

		if (this.session == null || id == null)
			return true;
		if (requestingSession.user == null)
			return requestingSession.equals(session);
		if (requestingSession.user.equals(session.user))
			return true;
		return (session.user == null || session.user.rank < requestingSession.user.rank);
	}

	@XmlElement(name="properties")
	public List<Property> getProperty()
	{
		// TODO: Call this only when required, prevent unnecessary marshalling cases
		if (Globals.getMarshallingOption(MarshallingOption.NO_BASKET_DETAILS))
			return null;

		if (modelProperties != null && !modelProperties.isEmpty() && Globals.session().contains(modelProperties.get(0)))
			return modelProperties;

		//TODO: Due to SERIOUS performance problems replaced fancy Criteria by native SQL query. Make it better someday... 
		String query = "select p.* from Basket b join BasketEntry be using (basket_id) join ExperimentalProperty ep using (exp_property_id) inner join Property p using (property_id) where b.basket_id = :id group by p.property_id";
		List<Property> props = Globals.session()
				.createSQLQuery(query)
				.addEntity(Property.class)
				.setLong("id", this.id)
				.list();

		// Group by parent
		List<Property> parents = new ArrayList<Property>();
		Iterator<Property> iter = props.iterator();
		while (iter.hasNext())
		{
			Property prop = iter.next();
			if (prop != null && prop.parent != null)
			{
				iter.remove();
				if (!parents.contains(prop.parent))
					parents.add(prop.parent);
			}
		}

		props.addAll(parents);

		return modelProperties = props;
	}

	@XmlAttribute(name="size")
	protected Long getBasketCountXML() 
	{
		if (Globals.getMarshallingOption(MarshallingOption.NO_BASKET_DETAILS)) // Significantly increase the load speed for large basket lists / Midnighter on Jun 17, 2011
			return null;
		if (cachedCount != null)
			return cachedCount;
		return getRowsSize();
	}

	@XmlTransient
	public Long getRowsSize() 
	{
		if (id != null && !Hibernate.isInitialized(entries))
		{
			logger.info("Counting records in basket " + name);
			return cachedCount = ((Long) Globals.session().createCriteria(BasketEntry.class).add(Restrictions.eq("basket", this)).setProjection(Projections.count("id")).uniqueResult()).longValue();
		}
		else 
			return (long) entries.size();
	}

	@XmlAttribute(name="models") // for XML only
	protected Long getModelsCount() 
	{
		if (!ThreadScope.get().controller.equals("basket") || Globals.getMarshallingOption(MarshallingOption.NO_BASKET_DETAILS))
			return null;
		return getModelCount(false);
	}

	@XmlAttribute(name="pending-models") // for XML only
	protected Long getPendingModelsCount() 
	{
		if (!ThreadScope.get().controller.equals("basket") || Globals.getMarshallingOption(MarshallingOption.NO_BASKET_DETAILS))
			return null;
		return getModelCount(true);
	}

	@XmlElementWrapper(name = "conditions-used")
	@XmlElement(name = "condition")
	protected List<Property> getConditions()
	{
		if (!Globals.getMarshallingOption(MarshallingOption.BASKET_CONDITIONS))
			return null;
		Criteria propertyCriteria = Globals.session().createCriteria(BasketEntry.class);
		propertyCriteria.add(Restrictions.eq("basket", this));
		propertyCriteria
		.createAlias("ep", "e")
		.createAlias("e.conditions", "c")
		.createAlias("c.values", "pv");

		return propertyCriteria.setProjection(Projections.groupProperty("pv.property")).list();
	}

	private Long getModelCount(boolean pending)
	{
		ModelsFilter filter = new ModelsFilter();
		filter.pendingTasks = pending;
		filter.trainingSets.add(this);

		return (Long) filter.createCriteria()
				.setProjection(Projections.countDistinct("id"))
				.uniqueResult();
	}

	public static void main(String[] args) {
		try
		{
			Globals.startAllTransactions();
			logger.info(Globals.session().createCriteria(Model.class)
					.createAlias("validationSets", "vs").
					add(Restrictions.eq("vs.id", 4106L)).setProjection(Projections.count("id"))
					.uniqueResult());
		}
		finally
		{
			Globals.rollbackAllTransactions();
		}
	}

	public Basket addEntries(Basket basket)
	{
		if (basket != null)
		{
			for (BasketEntry entry : basket.entries)
				this.entries.add(new BasketEntry(entry.ep));
			cachedCount = null;
		}
		return this;
	}

	/**
	 * Proper way to add new entries to the Basket
	 * @param ep
	 * @return
	 */

	public BasketEntry addEntry(ExperimentalProperty ep)
	{
		BasketEntry be = new BasketEntry();
		be.basket = this;
		be.ep = ep;
		entries.add(be);
		cachedCount = null;
		return be;
	}

	public String getFilenameCompatibleName()
	{
		String replacePattern = "[ \\-()=/<>.,\\\\+]";
		return name.replaceAll(replacePattern, "_");
	}

	public Basket()
	{
		if (lastModified == null)
			lastModified = new Timestamp(Calendar.getInstance().getTimeInMillis());	
	}

	//
	// Used for selection basket
	//

	public void addEntry(Long recordID)
	{
		Globals.session().createSQLQuery("insert ignore into BasketEntry (basket_id, exp_property_id) values (:basket_id, :ep_id)").setLong("basket_id", id).setLong("ep_id", recordID).executeUpdate();
		markModified();
	}

	@XmlElement(name = "excludedCount")
	public Long getExcludedCount()
	{
		if (!"basket".equals(ThreadScope.get().controller) || Globals.getMarshallingOption(MarshallingOption.NO_BASKET_DETAILS) || Globals.getMarshallingOption(MarshallingOption.BASKET_LIST_MODE))
			return null;
		return (Long) Globals.session().createCriteria(BasketEntry.class).add(Restrictions.eq("basket", this))
				.add(Restrictions.eq("exclude", Boolean.TRUE)).setProjection(Projections.count("id")).uniqueResult();
	}

	public void removeEntry(Long recordID)
	{
		Globals.session().createSQLQuery("delete ignore from BasketEntry where basket_id = :basket_id and exp_property_id = :ep_id").setLong("basket_id", id).setLong("ep_id", recordID).executeUpdate();
		markModified();
	}

	public void excludeOrIncludeEntry(Long recordID, boolean exclude)
	{
		Globals.session().createSQLQuery("update BasketEntry set exclude = :exclude where basket_id = :basket_id and exp_property_id = :ep_id")
		.setLong("basket_id", id)
		.setLong("ep_id", recordID)
		.setBoolean("exclude", exclude)
		.executeUpdate();
		markModified();
	}

	public long excludeOrIncludeEntries(List<Long> recordIDs, boolean exclude)
	{
		Long affectedCount = (Long) Globals.session().createSQLQuery("select count(*) c from BasketEntry where basket_id=:basketId and exp_property_id in (:list) and exclude != :exclude")
				.addScalar("c", LongType.INSTANCE)
				.setLong("basketId", id)
				.setParameterList("list", recordIDs)
				.setBoolean("exclude", exclude)
				.uniqueResult();
		Globals.session().createSQLQuery("update BasketEntry set exclude=:exclude where basket_id=:basketId and exp_property_id in (:list)")
		.setParameterList("list", recordIDs)
		.setLong("basketId", id)
		.setBoolean("exclude", exclude)
		.executeUpdate();

		markModified();

		return affectedCount;
	}

	public void includeAllEntries()
	{
		Globals.session().createSQLQuery("update BasketEntry set exclude=0 where basket_id=:basketId").setParameter("basketId", id).executeUpdate();
		markModified();
	}

	public void markModified()
	{
		cachedCount = null;
		if (id != null && id > 0)
			Globals.session().createSQLQuery("update Basket set last_modified = now() where basket_id=:basketId").setParameter("basketId", id).executeUpdate();
	}

	public void invertExclusion(List<Long> recordIDs)
	{
		Globals.session()
		.createSQLQuery("update BasketEntry set exclude = !exclude where basket_id=:basketId and exp_property_id in (:epIds)")
		.setParameter("basketId", id)
		.setParameterList("epIds", recordIDs)
		.executeUpdate();
		markModified();
	}

	public void delete()
	{
		if (id > 0)
		{
			// This is required optimisation, because otherwise deletion of records with > 10,000 records takes eternity...
			Globals.session().createSQLQuery("delete from BasketEntry where basket_id=:basketId").setLong("basketId", id).executeUpdate();
			entries.clear();
			Globals.session().delete(this);
		}
	}

	public List<Long> getRecordIDs()
	{
		return Globals.session().createSQLQuery("select exp_property_id from BasketEntry where basket_id=" + id + " order by exp_property_id").addScalar("exp_property_id", LongType.INSTANCE).list();
	}

	public void addEntries(List<Long> ids)
	{
		addEntries(ids, 0);
	}

	public void addEntries(List<Long> ids, int molOption)
	{
		SQLQuery query;
		switch (molOption)
		{
		case 1: 
			query = Globals.session().createSQLQuery("insert ignore into BasketEntry(basket_id, exp_property_id) select "+id+", min(ep.exp_property_id) from ExperimentalProperty ep join Molecule m using (molecule_id) where ep.exp_property_id in (:list) group by m.mapping2_id");
			break;
		case 2:
			query = Globals.session().createSQLQuery("insert ignore into BasketEntry(basket_id, exp_property_id) select "+id+", min(ep.exp_property_id) from ExperimentalProperty ep join Molecule m using (molecule_id) where ep.exp_property_id in (:list) group by m.mapping1_id");
			break;
		default:
			query = Globals.session().createSQLQuery("insert ignore into BasketEntry(basket_id, exp_property_id) select "+id+", exp_property_id from ExperimentalProperty where exp_property_id in (:list)");
			break;
		}
		query.setParameterList("list", ids).executeUpdate();
		markModified();
	}

	public void removeEntries()
	{
		Globals.session().createSQLQuery("delete ignore from BasketEntry where basket_id=" + id ).executeUpdate();
		markModified();
	}

	public void removeEntries(List<Long> ids)
	{
		Globals.session().createSQLQuery("delete ignore from BasketEntry where basket_id=" + id + " and exp_property_id in (:list)").setParameterList("list", ids).executeUpdate();
		markModified();
	}

	public boolean containsEntry(Long entryId)
	{
		List<Long> result = Globals.session().createSQLQuery("select count(*) cnt from BasketEntry where basket_id = :basket_id and exp_property_id = :ep_id").addScalar("cnt",LongType.INSTANCE).setLong("basket_id", id).setLong("ep_id", entryId).list();
		return (result.get(0) > 0);
	}

	public List<Long> getIds()
	{
		Criteria c = Globals.session().createCriteria(ExperimentalProperty.class)
				.createAlias("basketEntries", "bes")
				.add(Restrictions.eq("bes.basket", this))
				.setProjection(Projections.id());
		List<Long> results = c.list();
		return results;
	}


	@XmlElement // for XML only
	protected String getOwnerName()
	{
		if (user != null)
			return session == null || session.user == user? user.login : user.login +", created by "+ session.user.login;
		return null;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((id == null) ? 0 : id.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Basket other = (Basket) obj;
		if (id == null) {
			if (other.id != null)
				return false;
		} else if (!id.equals(other.id))
			return false;
		return true;
	}

	public void addMolecule(String sdf) throws IOException, TimeoutException
	{
		Molecule molecule = MoleculePeer.getMolecule(sdf);
		addMolecule(molecule);
	}


	public void addMolecule(Molecule mol)
	{
		ExperimentalProperty ep = new ExperimentalProperty();
		ep.molecule = mol;
		entries.add(new BasketEntry(ep));
	}

	private transient List<Long> ids;

	public List<Long> getIdentifiers()
	{
		if (ids == null)
		{
			ids = new ArrayList<Long>();
			for (BasketEntry entry : entries) {
				ids.add(entry.ep.id);
			}
		}

		return ids;
	}

	public static String getFreeBasketName(String name)
	{
		name = OCHEMUtils.getFilteredBasketName(name); // in case if previous name was using old standards

		int i = 0;
		String freeName;
		while (Basket.getBasket(freeName = name + (i == 0 ? "" : " " + i), Globals.userSession(), false) != null)
			i++;

		return freeName;
	}

	@Override
	public User getIntroducer() 
	{
		return user;
	}

	@Override
	public User getOwner() 
	{
		return user;
	}

	@Override
	public Integer getRights() 
	{
		return rights;
	}

	public List<PropertyOption> getPropertyOptions(Property pr) {
		Criteria propertyCriteria = Globals.session().createCriteria(BasketEntry.class);

		propertyCriteria.createAlias("ep", "e");
		propertyCriteria.add(Restrictions.eq("basket", this));
		propertyCriteria.add(Restrictions.eq("e.property", pr));

		propertyCriteria.setProjection(Projections.projectionList()
				.add(Projections.groupProperty("e.option")));

		List<PropertyOption> optionsList = propertyCriteria.list();

		return optionsList;
	}

	public void evict(boolean forced) {
		if (id != null || forced) // we do it only for "real baskets"; forced is required for baskets with records from the database
		{
			Globals.session().evict(this);
			for (BasketEntry be : entries) 
				if (Globals.session().contains(be))
					Globals.session().evict(be);
		}
	}

	private Set<Integer> stringToSet(String s) {
		Set<Integer> excludedBasketEntries = new HashSet<Integer>();
		String[] entries = s.split(",");
		for (String entry : entries)
			if (!entry.equals(""))
				excludedBasketEntries.add(Integer.valueOf(entry));
		return excludedBasketEntries;
	}


	private String setToString(Set<Integer> excludedBasketEntries) {
		if (excludedBasketEntries.size() == 0)
			return "";

		List<Integer> ids = new ArrayList<Integer>();
		ids.addAll(excludedBasketEntries);
		String s = ids.get(0).toString();

		for (int i=1; i<ids.size(); i++)
			s+=(","+ids.get(i));
		return s;
	}

	/**
	 * 
	 * @return Property_id mapping to mapping2_ids to be excluded 
	 */


	public Map<Long, Set<Integer>>  gexExcludedImplicitMolecules() {
		Map<Long,Set<Integer>> excludedImplicitMoleculesEntries = new HashMap<Long,Set<Integer>>();
		if(excludedrecords == null || excludedrecords.length() == 0) return excludedImplicitMoleculesEntries;
		for(String property: excludedrecords.split("\\R")) {
			String pieces[] = property.split(":");
			Long p = Long.parseLong(pieces[0]);
			Set<Integer> set = stringToSet(pieces[1]);
			excludedImplicitMoleculesEntries.put(p, set);
		}
		return excludedImplicitMoleculesEntries;
	}

	public String setExcludedImplicitMolecules(Map<Long, Set<Integer>> excludedImplicitMoleculesIds) {
		if (excludedImplicitMoleculesIds == null || excludedImplicitMoleculesIds.size() == 0)
			return excludedrecords = "";
		excludedrecords ="";
		for(Long property: excludedImplicitMoleculesIds.keySet()) {
			Set <Integer> propertySet = excludedImplicitMoleculesIds.get(property);
			if(!propertySet.isEmpty())
				excludedrecords += property + ":" + setToString(propertySet) +"\n";
		}
		return excludedrecords;		
	}

}

