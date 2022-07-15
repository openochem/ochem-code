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

package com.eadmet.business;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.ProjectionList;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Predicate;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.Tag;
import qspr.entities.Unit;
import qspr.entities.UnitCategory;
import qspr.entities.User;
import qspr.frontend.LabeledValue;
import qspr.workflow.utils.QSPRConstants
;
import qspr.util.AccessChecker;

import com.eadmet.business.PropertiesFilter.ApprovalStatus;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.useractions.EventFactory;
import com.eadmet.useractions.PropertyAction;
import com.eadmet.useractions.PropertyAction.PropertyActionType;
import com.eadmet.utils.mailer.Mailer;

@SuppressWarnings("unchecked")
public class PropertiesService 
{
	private ResultBuilder<Property> getResultBuilder(PropertiesFilter filter)
	{
		if (filter == null)
			filter = new PropertiesFilter();

		Criteria c = Globals.session().createCriteria(Property.class);

		if (filter.name != null)
			c.add(Restrictions.like("name", "%"+filter.name+"%"));
		else if (filter.query != null)
		{
			Disjunction d = Restrictions.disjunction();
			d.add(Restrictions.like("name", "%"+filter.query+"%"));
			d.add(Restrictions.like("aliases", "%"+filter.query+"%"));
			d.add(Restrictions.like("description", "%"+filter.query+"%"));
			c.add(d);
		}

		if (filter.id != null)
			c.add(Restrictions.eq("id", filter.id));
		else
			c.add(Restrictions.eq("isCondition", filter.condition));

		if (filter.parentId != null)
			c.createAlias("parent", "par").add(Restrictions.eq("par.id",filter.parentId));

		if (filter.directories != null)
			c.add(Restrictions.eq("isDirectory", filter.directories));

		if (filter.basketId != null)
		{
			Basket basket = (Basket) Globals.session().get(Basket.class, filter.basketId);
			c.createCriteria("experimentalProperties").createCriteria("basketEntries").add(Restrictions.eq("basket", basket));
		}

		if (filter.tagId != null)
			c.createAlias("tags", "t").add(Restrictions.eq("t.id", filter.tagId)); // Show properties from given tag
		else
			Globals.applyTaginationFilters(c, Property.class, "p"); // Apply global tagination filters

		if (filter.approvalStatus == ApprovalStatus.AWAITING_ONLY)
		{
			c.add(Restrictions.eq("approved", false));
			c.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE));
		}

		if (filter.name != null || filter.query != null)
			c.addOrder(Order.asc("name"));
		else
			c.addOrder(Order.desc("id"));

		AccessChecker.addAccessRestrictions(c, Globals.RIGHTS_FREELY_AVAILABLE, Globals.userSession().user, null,  filter.approvalStatus != ApprovalStatus.APPROVED_ONLY);

		return new ResultBuilder<Property>(c, Property.class, Globals.session()).setDistinct(true); //Properties are always distinct
	}


	public List<Property> get(PropertiesFilter filter, PaginationFilter pager)
	{
		ResultBuilder<Property> rb = getResultBuilder(filter).withPager(pager);
		pager.totalSize = rb.count();
		List<Property> properties = rb.list();
		return properties;
	}


	public List<Predicate> getPredicates()
	{
		return Globals.session().createCriteria(Predicate.class).list();
	}

	public List<UnitCategory> getUnitCategories()
	{
		return Globals.session().createCriteria(UnitCategory.class).list();
	}

	public List<LabeledValue> getLabels(PropertiesFilter filter, PaginationFilter pager)
	{
		List<Property> properties = get(filter, pager);
		List<LabeledValue> labels = new ArrayList<LabeledValue>();
		for (Property property : properties)
			labels.add(new LabeledValue(property.getName(), property.id));
		return labels;
	}

	public void delete(Property p)
	{
		AccessChecker.requestModificationPermission(p);

		// Delete unreferenced condition sets
		Globals.session().createSQLQuery("delete from PropertyValue where property_id=? and not exists (select * from ExperimentalProperty ep where ep.con_set_id = PropertyValue.con_set_id)")
		.setLong(0, p.id)
		.executeUpdate();

		Globals.session().createSQLQuery("delete from ArticleProperty where property_id=?")
		.setLong(0, p.id)
		.executeUpdate();

		Globals.session().delete(p);
		EventFactory.document("Property deletion", new PropertyAction(p, PropertyActionType.DELETE));
	}

	public Property approve(Property p)
	{
		AccessChecker.requestSuperuserPrivileges();
		p.approve();
		Globals.session().saveOrUpdate(p);
		return p;
	}

	public Property unapprove(Property p)
	{
		AccessChecker.requestSuperuserPrivileges();
		p.unapprove();
		Globals.session().saveOrUpdate(p);
		return p;
	}

	public Property addchild(Property parent, Property child)
	{
		if (!parent.isDirectory)
			throw new UserFriendlyException(parent.getName() + " is not a group");
		child.parent = parent;
		Globals.session().saveOrUpdate(child);
		return child;
	}

	public Property removechild(Property child)
	{
		AccessChecker.requestSuperuserPrivileges();
		child.parent = null;
		Globals.session().saveOrUpdate(child);
		return child;
	}

	public Property publish(Property p)
	{
		AccessChecker.requestModificationPermission(p);
		p.rights = Globals.RIGHTS_FREELY_AVAILABLE;
		p.updateHash();
		if (p.hasConflicts())
			throw new UserFriendlyException("Property with name \"" + p.getName() + "\" already exists! Please rename your property before publishing.");
		Globals.session().saveOrUpdate(p);
		Mailer.notifyAdmins("request to make public property " + p.getName(), " By user " + Globals.userSession().user.login);
		return p;
	}

	public Property applyTags(String propertyNameMask, String propertyTagName)
	{
		if (propertyNameMask != null && propertyTagName != null)
		{
			List<Property> properties = Globals.session().createCriteria(Property.class).add(Restrictions.like("name", "%"+propertyNameMask+"%")).list();
			for (Property property : properties) 
			{
				property.tags.add(Tag.getByName(propertyTagName));
				Globals.session().saveOrUpdate(property);
			}
		}
		return null;
	}

	public Property create(PropertiesAction action)
	{
		if (Globals.isGuestUser())
			throw new UserFriendlyException("Guest users cannot create properties. Please, register (its free) to access extended functionality.");

		Property p = new Property();
		p.introducer = Globals.userSession().user;
		p.isCondition = action.isCondition;
		p.isDirectory = action.isDirectory;
		if (p.isCondition)
		{
			p.publish();
			p.approved = true;
		}
		return p;
	}


	public Property saveOptions(Property property, String[] optionIds, String[] optionNames)
	{
		property.owner = Globals.userSession().user;

		if (optionIds == null)
		{
			property.options.clear();
			Globals.session().saveOrUpdate(property);
			return property;
		}

		List<PropertyOption> savedOptions = new ArrayList<PropertyOption>();
		List<String> savedNames = new ArrayList<String>();

		for (int i = 0; i < optionIds.length; i++)
		{
			Long id = Long.valueOf(optionIds[i]);
			PropertyOption co = null;

			if (id > 0) 
				co = (PropertyOption) Globals.session().get(PropertyOption.class, id);
			else {

				if(!PropertyOption.validateOption(optionNames[i]))
					throw new UserFriendlyException("Property options should be a non-empty String and not a Float/Integer, but you provided: \"" + optionNames[i] + "\"");

				co = Repository.option.getPropertyOptionByName(optionNames[i], property.id, true, false);
			}

			/*// To check whether it works correctly ... id <= 0 corresponds to the above conditions
			if (id < 0 || !co.isReferenced())
			{
				// POLICY: Only unreferenced options are editable
				co.property = property;
				co.name = optionNames[i].trim();
			}
			 */
			if (!savedNames.contains(co.name))
			{
				savedOptions.add(co);
				savedNames.add(co.name);
			}
		}
		property.options.clear();
		property.options.addAll(savedOptions);
		Globals.session().saveOrUpdate(property);
		return property;
	}

	public Property edit(PropertiesAction action)
	{
		Property p = null;

		if (action.id == null || (action.id == -1))
			p = create(action);
		else
			p = (Property)Globals.session().get(Property.class, action.id);

		AccessChecker.requestModificationPermission(p);

		p.owner = Globals.userSession().user;

		if (p.introducer == null)
			p.introducer = Globals.userSession().user;

		if (p.getPropertyCount() == 0)
			if (action.name != null)
				p.setName(action.name);

		if (p.id == null)
		{
			p.type = action.propertyType;

			if (p.isNumeric())
				p.unitCategory = (UnitCategory) Globals.session().get(UnitCategory.class, action.unitCategoryId);
			else
			{
				p.unitCategory = UnitCategory.getByName(QSPRConstants.CLASS);
				p.defaultUnit = p.unitCategory.getDefaultUnit();
			}

			if (p.isTextual() && !p.isCondition)
				throw new UserFriendlyException("Arbitrary text properties are not currently supported");
		}
		else
		{
			UnitCategory newCategory = (UnitCategory) Globals.session().get(UnitCategory.class, action.unitCategoryId);
			if (!newCategory.equals(p.unitCategory))
				updateUnitCategory(p, newCategory, action.confirmed);
		}

		if (action.bonusPointsWeight != null && Globals.isSuperUser())
			p.bonusPointsWeight = action.bonusPointsWeight;

		if (action.parentId != null)
			p.parent = (Property) Globals.session().get(Property.class, action.parentId);

		p.description = action.description;
		p.aliases = action.aliases;

		if (p.isNumeric())
		{
			if (action.defaultUnitId != null)
				p.defaultUnit = (Unit) Globals.session().get(Unit.class, action.defaultUnitId);

			if (p.unitCategory == null)
				throw new UserFriendlyException("Unit category can not be empty for numeric properties");
			if (p.defaultUnit == null)
				throw new UserFriendlyException("Default unit can not be empty for numeric properties");
		}

		if (action.isPublic)
		{
			if (Globals.userSession().user == null)
				throw new UserFriendlyException("Guest users cannot create public records. Please, register and login as a registered user to introduce public data.");

			p.publish();
		}

		if (action.isApproved && !p.approved)
		{
			AccessChecker.requestModeratorPrivileges();
			p.approve();
		}

		Set<Long> newConditionIDs = new HashSet<Long>();
		newConditionIDs.addAll(action.obligatoryConditionIds);

		Set<Long> oldConditionIDs = new HashSet<Long>();
		for (Property condition : p.obligatoryConditions)
			oldConditionIDs.add(condition.id);

		if (!oldConditionIDs.equals(newConditionIDs))
		{
			long recordReferences = 0;
			if (!p.isCondition && !p.isDirectory && (p.id != null))
			{
				Criteria c = Globals.session().createCriteria(ExperimentalProperty.class)
						.add(Restrictions.eq("property", p))
						.setProjection(Projections.countDistinct("id"));
				recordReferences = (Long)c.list().get(0);
			}

			if (recordReferences == 0 || (Globals.userSession().user != null && (Globals.userSession().user.equals(p.introducer) || Globals.userSession().user.isSuperUser())))
			{
				int changes = Globals.session().createSQLQuery("update ExperimentalProperty set ep_md5 = null where ((deleted is null) and (property_id = "+p.id+"))").executeUpdate();
				logger.info("[Obligatory Conditions Change] Update of obligatory conditions invalidated hashes for "+changes+" records. They hopefully will get recalculated by the cronjob.");
				changes = Globals.session().createSQLQuery("update ExperimentalProperty set status = null where ((deleted is null) and (status is not null) and (property_id = "+p.id+"))").executeUpdate();
				logger.info("[Obligatory Conditions Change] Additionally status was changed from invalid to valid for "+changes+" records. They may stay valid after hash recalculation.");
				p.obligatoryConditions.clear();

				for (Long conditionId : action.obligatoryConditionIds) 
					p.obligatoryConditions.add((Property)Globals.session().get(Property.class, conditionId));
			}
			else
				throw new UserFriendlyException("You are not authorized to update obligatory condition for this property, since it already has data associated with it");
		}

		p.updateHash();
		if (p.hasConflicts())
			throw new UserFriendlyException("Property with name \"" + p.getName() + "\" already exists! Please, use this property or choose another name");

		if (p.id == null)
			EventFactory.document("Property creation", new PropertyAction(p, PropertyActionType.CREATE));

		Globals.session().saveOrUpdate(p);
		return p;
	}

	private void updateUnitCategory(Property property, UnitCategory newCategory, boolean confirmed)
	{
		// A set of actions should be done to change the system of units of a property
		// Every unit, associated with this property should be:
		//   eithier (a) moved to the new system of units, if no other properties refer to this unit
		//   or (b) be copied to the new system of units, preserving the unit name
		// This is required for orderinng of the EMBL data
		// Midnighter

		// First check if there is hidden data with this property. If yes, the system of units cannot be changed
		List<User> introducers = Globals.session().createCriteria(ExperimentalProperty.class)
				.add(Restrictions.eq("property", property))
				.add(Restrictions.ne("rights", Integer.valueOf(Globals.RIGHTS_FREELY_AVAILABLE)))
				.add(Restrictions.ne("introducer", Globals.userSession().user))
				.setProjection(Projections.groupProperty("introducer")).list();
		if (introducers.size() > 0)
			throw new UserFriendlyException("Unit category cannot be changed, since there is hidden data for this property, introduced by " + introducers);


		StringWriter stringWriter = new StringWriter();
		PrintWriter writer = new PrintWriter(stringWriter);

		writer.println("Dear user,\n\n you are going to change the system of units for the property \"" + property + "\" from \"" + property.unitCategory.name + "\" to \"" + newCategory.name + "\"");
		writer.println("This will require following actions:\n");

		// Find all the units, used with this property
		List<Unit> units = Globals.session()
				.createCriteria(ExperimentalProperty.class).add(Restrictions.eq("property", property))
				.add(Restrictions.isNotNull("unit"))
				.setProjection(Projections.groupProperty("unit"))
				.list();

		// Iterate over all the detected units and decide whether we need to dublicate it or just move to the new unit category
		for (Unit unit : units) 
		{
			// What are the other referencing properties?
			List<Property> properties = Globals.session().createCriteria(ExperimentalProperty.class)
					.add(Restrictions.eq("unit", unit))
					.add(Restrictions.ne("property", property))
					.setProjection(Projections.groupProperty("property"))
					.list();

			writer.println(unit.getName() + ": " + properties.size() + " other referencing properties");
			Unit newUnit = newCategory.getUnitByName(unit.getName());

			if (newUnit == null && properties.size() == 0)
			{
				writer.println("** Moving \"" + unit + "\" to \"" + newCategory.name + "\"");
				if (confirmed)
				{
					unit.category = newCategory;
					Globals.session().saveOrUpdate(unit);
				}
			}
			else
			{
				if (newUnit == null)
				{
					writer.println("** Creating new unit \""+unit.getName()+"\" in category \"" + newCategory + "\"");
					newUnit = new Unit();
					newUnit.category = newCategory;
					newUnit.introducer = newUnit.owner = Globals.userSession().user;
					newUnit.setName(unit.getName());
					if (confirmed)
						Globals.session().saveOrUpdate(newUnit);
				}


				Long count = (Long) Globals.session().createCriteria(ExperimentalProperty.class)
						.add(Restrictions.eq("unit", unit))
						.add(Restrictions.eq("property", property))
						.setProjection(Projections.count("id")).uniqueResult();


				writer.println("** Changing unit for " + count +" records to \"" + unit + "\" (" + newCategory + ")");
				if (confirmed)
				{
					// In this case HQL query seems the only option to update records fast
					Globals.session().createQuery("update ExperimentalProperty set unit=:newUnit where unit=:oldUnit and property=:property")
					.setParameter("oldUnit", unit)
					.setParameter("newUnit", newUnit)
					.setParameter("property", property)
					.executeUpdate();
				}


				if (properties.size() == 0)
				{
					writer.println("** Deleting unit \"" + unit + "\" from \""+unit.category+"\"");
					if (confirmed)
						Globals.session().delete(unit);
				}
			}
		}

		// At last, after all the stuff, change the system of units
		property.unitCategory = newCategory;

		if (!confirmed)
			throw new UserFriendlyException(stringWriter.toString());
		else
			logger.info(stringWriter.toString());
	}

	@SuppressWarnings("rawtypes")
	public void fillUsedConditions(Property prop)
	{
		if (prop.id != null && !prop.isCondition)
		{
			ProjectionList projList = Projections.projectionList();
			projList.add(Projections.groupProperty("val.property"));
			projList.add(Projections.alias(Projections.countDistinct("id"),"cnt"));

			Criteria c = Globals.session().createCriteria(ExperimentalProperty.class)
					.add(Restrictions.isNull("deleted"))
					.createAlias("conditions", "con")
					.createAlias("con.values", "val")
					.add(Restrictions.eq("property",prop))
					.setProjection(projList)
					.addOrder(Order.desc("cnt"));

			ExperimentalProperty.addAccessRestrictions(c, Globals.RIGHTS_FREELY_AVAILABLE, Globals.userSession().user, false, true);

			List conditionsUsed = c.list();

			if (conditionsUsed.size() > 0)
			{
				prop.conditionsUsed = new ArrayList<Property>();
				for (int i=0; i<conditionsUsed.size(); i++)
				{
					Object[] tuple = (Object[])conditionsUsed.get(i);
					Property condition = (Property)tuple[0];
					condition.count = (Long)tuple[1];
					prop.conditionsUsed.add(condition);
				}
			}
		}
	}


	private static Logger logger = LogManager.getLogger(PropertiesService.class);
}
