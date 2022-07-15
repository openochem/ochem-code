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

package qspr.controllers;

import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.Query;
import org.hibernate.criterion.Conjunction;
import org.hibernate.criterion.DetachedCriteria;
import org.hibernate.criterion.Disjunction;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.criterion.Subqueries;
import org.hibernate.sql.JoinType;
import org.hibernate.type.IntegerType;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.ThreadScope;
import qspr.business.DataApprovalEmailNotifier;
import qspr.business.SubstructureSearchFilter;
import qspr.business.WebFilters;
import qspr.dao.Repository;
import qspr.entities.Alert;
import qspr.entities.Article;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Mapping1;
import qspr.entities.Mapping2;
import qspr.entities.Mapping2Filter;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.Molecule;
import qspr.entities.Predicate;
import qspr.entities.Property;
import qspr.entities.Session;
import qspr.entities.Tag;
import qspr.entities.User;
import qspr.exception.DublicateRestoreException;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.modelling.ModelStatistics;
import qspr.modelling.PointSelector;
import qspr.modelling.PointStatistics;
import qspr.modelling.SetStatistics;
import qspr.toxicity.RequestStatistics;
import qspr.util.AccessChecker;
import qspr.util.CASRN;
import qspr.util.CriteriaWrapper;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.business.ModelDotService;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.messaging.DialogueService;
import com.eadmet.useractions.BasketAction;
import com.eadmet.useractions.BasketAction.BasketActionType;
import com.eadmet.useractions.EventFactory;
import com.eadmet.utils.OCHEMUtils;

/**
 * The controller for the experimental property browser.
 * Requires significant refactoring; first of all - separation of the filtering logic to a separate class similarly to AlertsFilter.java or ECommerceMoleculesFilter.java
 * 
 * @author all, really everyone
 *
 */
@SuppressWarnings("unchecked")
@Controller
public class EPBrowserController extends BrowserWrapper
{
	private static final String ACTION_PUBLISH_FREE = "publish_freely_available";

	public Map<String, Long> md5Hashes = new HashMap<String, Long>();
	private Map<Long, List<Long>> recordsDeletedFromModel = new HashMap<Long, List<Long>>();
	private Map<Long, Integer> basketToMapping2FilterCache = new HashMap<Long, Integer>(); 

	@Autowired
	public DialogueService dialogueService;

	public EPBrowserController()
	{
		sessionRequired = true;
	}

	private Criteria createConditionCriteria(WebFilters filters, int rights, Criteria criteria)
	{
		List<Long> conditionIds = new ArrayList<Long>();
		Map<String, Object> parameters = new HashMap<String, Object>();

		for (String filter : filters.getFilterSet())
			if (filter.startsWith("cond-id"))
				conditionIds.add(Long.valueOf(filters.get(filter)));				

		StringBuilder query = new StringBuilder("select distinct cs from ConditionSet cs");
		if (conditionIds.size() > 0)
		{
			for (Long key : conditionIds)
				query.append(" left join cs.values as vals"+key);

			query.append(" where (1=1)");

			for (Long key : conditionIds)
			{
				query.append(" and ((vals"+key+".property.id = "+key+")");//3
				Property condition = (Property) Globals.session().get(Property.class, key);

				if (condition.isQualitative())
					if (filters.notEmpty("cond-option" + key))
						query.append(" and (vals"+key+".option.id = " + filters.getLong("cond-option" + key) + ")");

				if (condition.isNumeric())
					if (filters.notEmpty("cond-value" + key))
					{
						query.append(" and (vals"+key+".unit.id = " + filters.getLong("cond-unit" + key) + ")");//1
						query.append(" and (vals"+key+".value = :value" + key + ")");
						parameters.put("value" + key, Double.valueOf(filters.get("cond-value" + key)));
					}

				if (condition.isTextual())
					if (filters.notEmpty("cond-text-value" + key))
					{
						query.append(" and (vals"+key+".textualValue like :tvalue" + key + ")");
						parameters.put("tvalue" + key, filters.get("cond-text-value" + key) + "%");
					}

				query.append(")");
			}
			Query ourHql = Globals.session().createQuery(query.toString());

			for (String key : parameters.keySet())
				ourHql.setParameter(key, parameters.get(key));

			List<Object> tmp = ourHql.list();
			if (tmp.size() > 0)
				criteria.add(Restrictions.in("conditions", tmp));
			else
				criteria.add(Restrictions.sqlRestriction("1 = 0"));
		}



		return criteria;
	}

	private List<Long> getSelectedIds(PointSelector pointSelector, ModelStatistics ms, Integer set, Long real, Long predicted){

		SetStatistics ss = ms.sets.get(set);

		if(pointSelector!= null && !ss.setId.equals(QSPRConstants.TRAINING)){ // we have a problem and need to init pointSelector first
			SetStatistics sss = ms.sets.get(0);
			pointSelector.pointMatches(sss.points.get(0), sss);
		}

		List<Long> ids = new ArrayList<Long>();
		for (PointStatistics ps : ss.points) {
			if (Math.round(ps.real) == real && Math.round(ps.predicted) == predicted && ps.error == null)
				if (pointSelector == null || pointSelector.pointMatches(ps, ss))
					ids.add(Math.abs(ps.id));
		}
		return ids;
	}


	private void excludeRecords(Property prop, Basket basket, List<Long> ids) {
		{
			if(ids.size() == 0) return;

			Map<Long,Set<Integer>> entries = basket.gexExcludedImplicitMolecules();
			HashSet <Integer> set = new HashSet<Integer>();
			entries.put(prop.id, set);
			for(Long id: ids) 
				set.add(Repository.record.getRecord(id).molecule.mapping2.id);
			basket.setExcludedImplicitMolecules(entries);

			Globals.session().saveOrUpdate(basket);

			//throw new UserFriendlyException(basket.description);
		}

	}


	protected CriteriaWrapper createMainCriteria(WebFilters filters, int rights, Criteria criteria, Mutable<Boolean> distinctRequired) throws Exception
	{
		boolean distinctQueryRequired = false;
		User user = Globals.userSession().user;
		//Condition filters parse block
		criteria = createConditionCriteria(filters, rights, criteria);
		CriteriaWrapper wrapper = new CriteriaWrapper(criteria);

		// How to get filters.has("id") to true?!?
		if (filters.has("id"))
		{
			// ID(s) have been explicitly specified
			// No need to look for other filters
			filterById(criteria);
		}
		else
		{		
			if (filters.has("tag"))
			{
				Set<Tag> tagFilters = Globals.getTaginationFilters(null);
				Tag tag = (Tag)Globals.session().get(Tag.class, getLongParam("tag"));
				if (!tagFilters.contains(tag))
					tagFilters.add(tag);
			}

			// Some actions always for on a selection!
			if (filters.in("action", new String[]{"addbasket", "removebasket", "deleteselected", "addtag", "removetag", "removeselect", "restoreselected", 
					"approve", "disapprove", "unapprove", "publish", ACTION_PUBLISH_FREE}))
			{
				//Set<Long> selectionList = Globals.userSession().selectionList;
				Basket selectionBasket = Globals.selectionBasket(assertParam("trash"));
				if (selectionBasket.getRowsSize() == 0) // including "-1"
					throw new UserFriendlyException("You did not select any records");
				else
				{
					wrapper.createAlias("selectedBasketEntries","selBasketEntries");
					criteria.add(Restrictions.eq("selBasketEntries.basket", selectionBasket));
					//criteria.add(Restrictions.in("id", selectionList));
				}
			}

			if (filters.has("id1") || filters.has("name") || filters.has("duplicates") ||  filters.has("mmlo") || filters.has("mmhi")
					|| filters.has("inchierror") || filters.has("emptymol") 
					|| filters.has("compare-baskets") || filters.has("uniquemol") || filters.has("cadasterSearch")
					|| (
							filters.has("sortby") 
							&& (filters.get("sortby").equals("molecule") || filters.get("sortby").equals("molprop") || filters.get("sortby").equals("mol.molWeight"))
							)
					|| Globals.getTaginationFilters(Mapping2.class).size() > 0
					|| filters.has("xemistry-smiles"))
			{
				wrapper.createAlias("molecule", "mol");
				wrapper.createAlias("mol.mapping1", "mp1");
				wrapper.createAlias("mol.mapping2", "mp2");
			}

			SubstructureSearchFilter ssFilter = new SubstructureSearchFilter(filters);
			boolean globalSSSearch = false;
			if (ssFilter.structure != null)
			{
				// Some heuristics: there are two variants for substructure search
				if (filters.has("basket-select") && (filters.getInteger("basket-select") > 0) || filters.has("property") || filters.has("article"))
				{
					// Variant 1: If property or a basket is selected, use Xemistry UDF
					wrapper.createAlias("mp2.xemistryIndex", "xemistryIndex");
					if (ssFilter.isSubstructureSearch())
						criteria.add(Restrictions.sqlRestriction("match_substructure(screen, molecule, '" + ssFilter.structure + "')"));
					else
					{
						String query = "similarity('" + ssFilter.structure + "', simscreen) > " + ssFilter.similarityThreshold;
						criteria.add(Restrictions.sqlRestriction(query));
						// TODO: Consider ordering the results by similarity
						//order = " order by  similarity('" + smiles + "',screen) desc";
					}
				}
				else
				{
					// Variant 2: If its a global search, we save the matching molecules in a Mapping2 filter

					// Do not allow this at the moment. Its a too unspecific search.
					throw new UserFriendlyException("To use substructure search, please make your search more specific by activating one of the following filters: property, article or basket");

					//globalSSSearch = true;
					//wrapper.createAlias("mp2.filters", "mp2Filters");
					//criteria.add(Restrictions.eq("mp2Filters.filterId", ssFilter.getFilterID()));
				}
			}

			if (filters.has("cs"))
				criteria.add(Restrictions.eq("conditions.id", filters.getLong("cs")));

			if (filters.has("mmlo"))
				criteria.add(Restrictions.ge("mol.molWeight", new Double(filters.get("mmlo"))));

			if (filters.has("mmhi"))
				criteria.add(Restrictions.le("mol.molWeight", new Double(filters.get("mmhi"))));

			if (filters.has("property") || filters.has("dummy") || filters.has("hidedummy")
					|| filters.has("duplicates") || filters.has("uniquemol")
					|| filters.has("compare-baskets") 
					|| Globals.getTaginationFilters(Property.class).size() > 0)
				wrapper.createAlias("property", "p");

			if (filters.has("predicate"))
				criteria.add(Restrictions.eq("predicate", Predicate.get("=")));

			if (filters.has("condition"))
			{
				Property condition = (Property) Globals.session().get(Property.class, getLongParam("condition"));
				criteria.createCriteria("conditions").createCriteria("values").add(Restrictions.eq("property", condition));
				distinctQueryRequired = true;
			}

			if (filters.has("structure-basket") && filters.has("structure-basket-dofilter") && (filters.getLong("structure-basket") != -1) && (filters.getLong("structure-basket-dofilter") != 0))
			{
				/*
				 * We reuse the Mapping2Filter functionality to join molecules from a basket.
				 * We keep a basket->filter map, but only to clean up the old filter on repetitive queries.
				 * Reusing the filter for repetitive queries is theoretically possible, but somewhat difficult,
				 * as we would need to watch closely when the basket has changed and invalidate the filter.
				 * Not necessary right now, as deleting/creating filters on the fly is relatively fast even for huge baskets.
				 * 
				 * NoS 05.02.2014
				 */
				long timer;
				Basket b = (Basket) Globals.session().get(Basket.class, getLongParam("structure-basket"));
				Integer filterId = basketToMapping2FilterCache.get(b.id);

				timer = System.nanoTime();
				if (filterId != null) 
					Mapping2Filter.clear(filterId);
				logger.info("Cleared previous filter in "+(System.nanoTime() - timer)/1000000+"ms");

				timer = System.nanoTime();
				filterId = Mapping2Filter.generateFilterID();
				basketToMapping2FilterCache.put(b.id, filterId);
				Mapping2Filter.addBasketCompoundsToFilter(filterId, b.id);
				logger.info("Created new filter in "+(System.nanoTime() - timer)/1000000+"ms");

				wrapper.createAlias("molecule", "mol");
				wrapper.createAlias("mol.mapping1", "mp1");
				wrapper.createAlias("mp1.mapping2", "mp2");
				wrapper.createAlias("mp2.basketFilters", "bmp2Filters");
				criteria.add(Restrictions.eq("bmp2Filters.filterId", filterId));
				distinctQueryRequired = true;
			}

			if (filters.has("unit") && !filters.get("unit").equals("-1"))
			{
				criteria.createAlias("conditions", "c", JoinType.LEFT_OUTER_JOIN)
				.createAlias("c.values", "cv", JoinType.LEFT_OUTER_JOIN)
				.createAlias("cv.unit", "cu", JoinType.LEFT_OUTER_JOIN)
				.createAlias("cv.property", "cp", JoinType.LEFT_OUTER_JOIN);

				Conjunction unitInCondition = Restrictions.conjunction();
				unitInCondition.add(Restrictions.eq("cu.id", filters.getLong("unit")));
				unitInCondition.add(Restrictions.eq("cp.type", Property.TYPE_NUMERIC));
				criteria.add(Restrictions.or(Restrictions.eq("unit.id", filters.getLong("unit")), unitInCondition));
				distinctQueryRequired = true;
			}

			if (filters.has("mm_id") && filters.has("num") && filters.has("real") && filters.has("predicted") && filters.has("set"))
			{
				// TODO: Optimise in future - store the list of IDs in session, not to recalculate it every time
				ModelMapping modelMappings = (ModelMapping) Globals.session().get(ModelMapping.class, filters.getLong("mm_id"));
				ModelStatistics ms = null;
				if(filters.getInteger("num") == 1)
					ms	= (ModelStatistics) modelMappings.statisticsOriginal.getObject();
				else
					ms	= (ModelStatistics) modelMappings.statisticsRecalculated.getObject();

				if (assertParam("validationSetNum"))
					ms.validationSetId = getParam("validationSetNum");

				ms.actualizeStatistics(modelMappings);

				PointSelector pointSelector = ModelController.getPointSelector(ThreadScope.get().getHttpServletRequest());

				int set = filters.getInteger("set");

				List<Long> ids = getSelectedIds(pointSelector, ms, set, filters.getLong("real"), filters.getLong("predicted"));

				criteria.add(Restrictions.in("id", ids));

				if (ids.size() > 0 && set == 0 && ThreadScope.get().getHttpServletRequest().getParameter("exclude-records") != null) 
					excludeRecords(modelMappings.property,Repository.basket.getById(ms.sets.get(0).basketId), ids);

			}

			if (filters.has("mol-filter"))
			{
				wrapper.createAlias("molecule", "mol");
				wrapper.createAlias("mol.mapping2", "mp2");
				wrapper.createAlias("mp2.filters", "mp2Filters");
				criteria.add(Restrictions.eq("mp2Filters.filterId", getIntParam("mol-filter")));
			}

			if (filters.has("bu"))
			{
				criteria.createAlias("batchUploadRows", "bur");
				criteria.createAlias("bur.batchUpload", "bu");
				criteria.add(Restrictions.eq("bu.id", filters.getLong("bu")));
				distinctQueryRequired = true;
			}

			if (filters.has("article"))
				criteria.add(Restrictions.eq("article.id", filters.getLong("article")));
			if (filters.has("property"))
			{
				if (filters.get("property").contains(","))
				{
					String[] propertyStr = filters.get("property").split(",");
					Long[] propertyIDs = new Long[propertyStr.length];
					for (int i = 0; i < propertyIDs.length; i++)
						propertyIDs[i] = Long.valueOf(propertyStr[i]);
					criteria.add(Restrictions.in("p.id", propertyIDs));
				}
				else
				{
					Property property = (Property) Globals.session().get(Property.class, filters.getLong("property"));
					if (property.isDirectory)
						criteria.add(Restrictions.eq("p.parent", property));
					else
						criteria.add(Restrictions.eq("p.id", filters.getLong("property")));
				}
			}

			if (filters.has("hidedummy"))
			{
				Property property = Property.getByName(QSPRConstants.DUMMY);
				criteria.add(Restrictions.ne("p.id", property.id));
			}

			if (filters.has("dummy"))
			{
				Property property = Property.getByName(QSPRConstants.DUMMY);
				criteria.add(Restrictions.eq("p.id", property.id));
				filters.addFilter("property", property.id.toString(), QSPRConstants.DUMMY);
			}

			if (filters.has("option"))
				criteria.add(Restrictions.eq("option.id", filters.getLong("option")));

			if (filters.has("idarticle"))
				criteria.add(Restrictions.eq("artMolId", filters.get("idarticle")));
			if (filters.has("page"))
				criteria.add(Restrictions.eq("artPageNum", filters.getInteger("page")));
			if (filters.has("line"))
				criteria.add(Restrictions.eq("artLineNum", filters.getInteger("line")));			
			if (filters.has("table"))
				criteria.add(Restrictions.eq("artTableNum", filters.get("table")));	
			if (filters.has("id1"))
				criteria
				.add(Restrictions.eq("mp1.id", filters.getLong("id1")));

			if (filters.has("recordsDeletedFromModel"))
			{
				ModelMapping mm = (ModelMapping) Globals.session().get(ModelMapping.class, filters.getLong("recordsDeletedFromModel"));

				if (recordsDeletedFromModel.get(mm.id) != null)
					criteria.add(Restrictions.in("id", recordsDeletedFromModel.get(mm.id)));
				else
				{
					ModelStatistics ms = (ModelStatistics) mm.statisticsOriginal.getObject();
					ms.actualizeStatistics(mm);
					criteria.add(Restrictions.in("id", ms.deletedRecordsIds));
					recordsDeletedFromModel.put(mm.id, ms.deletedRecordsIds);	
				}
			}

			if (filters.has("similarto"))
			{
				Molecule similarMol = (Molecule) 
						Repository.molecule.getMolecule(filters.getLong("similarto"));
				criteria.add(Restrictions.eq("molecule.mapping1", similarMol.mapping1));			
			}

			if (filters.has("private"))
				criteria.add(Restrictions.eq("rights", Integer.valueOf(0)));

			if(filters.has("interval"))
			{
				if(filters.get("interval").contains(","))
				{
					String [] rec_time = filters.get("interval").trim().split(",");
					ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, Long.valueOf(rec_time[0].trim()));
					if(ep != null)
					{
						Long time = Long.valueOf(rec_time[1].trim());
						Timestamp timeBefore = new Timestamp(ep.time.getTime() - time*60L*1000L);
						Timestamp timeAfter = new Timestamp(ep.time.getTime() + time*60L*1000L);
						criteria.add(Restrictions.between("time", timeBefore, timeAfter));
					}
				}else
				{
					Timestamp timeNow = new Timestamp(Calendar.getInstance().getTimeInMillis()); 
					Timestamp timeBefore = new Timestamp(Calendar.getInstance().getTimeInMillis() - filters.getLong("interval")*60L*1000L);
					criteria.add(Restrictions.between("time", timeBefore, timeNow));
				}
			}
			if (filters.has("experimental"))
				criteria.add(Restrictions.eqProperty("connectedProperty.id", "id"));
			if (filters.has("validated"))
				criteria.add(Restrictions.eq("ep_status", ExperimentalProperty.STATUS_TOVERIFY));

			if (filters.has("inchierror"))
			{
				criteria.add(Restrictions.like("mp1.inchi1", "________________________________"));
				criteria.add(Restrictions.ne("mp1.inchi1", QSPRConstants.EMPTY_MD5));
			}

			if (filters.has("nameerror"))
			{
				distinctQueryRequired = true;
				if (filters.has("yesstereo"))
					criteria.createCriteria("colorednames").add(Restrictions.or(Restrictions.eq(QSPRConstants.VALIDATION, Long.valueOf("2")), Restrictions.eq(QSPRConstants.VALIDATION, Long.valueOf("4"))));
				else 
					criteria.createCriteria("colorednames").add(Restrictions.eq(QSPRConstants.VALIDATION, Long.valueOf("4")));
			}

			if (filters.has("emptymol"))
				criteria.add(Restrictions.eq("mp1.inchi1", QSPRConstants.EMPTY_MD5));

			if (filters.has("error"))
				criteria.add(
						Restrictions.in("ep_status", new Integer[]{ExperimentalProperty.STATUS_ERROR, ExperimentalProperty.STATUS_INVALID}));

			if (filters.has("primaryrecords"))
				criteria.add(Restrictions.eqProperty("firstEntry", "id"));
			if (filters.has("secondaryrecords"))
				criteria.add(Restrictions.eq("isPrimary", false));


			if (filters.has("introducer"))
			{
				if (filters.get("introducer").equals("group"))
				{
					if (user != null && user.group != null)
						criteria.add(Restrictions.eq("intr.group", user.group));
				}
				else if (filters.get("introducer").equals("modified-by-me"))
				{
					// Records that were modified, but not introduced by me
					User me = Globals.userSession().user;
					criteria.add(Restrictions.eq("owner", me));
					criteria.add(Restrictions.ne("introducer", me));
				}
				else if (filters.get("introducer").equals("modified-by-others"))
				{
					// Records that were introduced by me, but modified by others
					User me = Globals.userSession().user;
					criteria.add(Restrictions.eq("introducer", me));
					criteria.add(Restrictions.ne("owner", me));
				}
				else if (filters.get("introducer").equals("my-records"))
				{
					// Records that were introduced or modified by me
					User me = Globals.userSession().user;
					Disjunction authorRestriction = Restrictions.disjunction();
					authorRestriction.add(Restrictions.eq("owner", me));
					authorRestriction.add(Restrictions.eq("introducer", me));
					criteria.add(authorRestriction);
				}
				else
				{
					// Records, that were introduced by a specified user
					User filteredUser = filters.get("introducer").matches("[0-9]+") ? (User) Globals.session().get(User.getCurrentClass(), filters.getLong("introducer")) : User.getByLogin(filters.get("introducer"));
					//Disjunction authorRestriction = Restrictions.disjunction();
					//authorRestriction.add(Restrictions.eq("owner", filteredUser));
					criteria.add(Restrictions.eq("introducer", filteredUser));
					//criteria.add(authorRestriction);
				}

			}

			if (filters.has("hidden"))
				criteria.add(Restrictions.eq("rights", Globals.RIGHTS_NONE));
			if (filters.has("value"))
				criteria.add(Restrictions.eq("value", new Double(filters.get("value"))));
			if (filters.has("canonicalValue"))
				criteria.add(Restrictions.eq("canonicalValue", Double.valueOf(filters.get("canonicalValue"))));

			Globals.applyTaginationFilters(criteria, Property.class, "p");
			Globals.applyTaginationFilters(criteria, Mapping1.class, "mp1");

			// Molecule name filter
			if (filters.has("name"))
			{
				// Check if its number -> the search by QID (ID2)
				// Check if its inchie-key ->then..
				// Otherwise search by name
				String enteredName = filters.get("name");
				if (enteredName.matches("(M[0-9]+,?)+")) // MXXXXX (mapping2)
				{
					if (!enteredName.contains(","))
						criteria.add(Restrictions.eq("mp2.id", Integer.valueOf(enteredName.substring(1))));
					else
					{
						String[] vals = enteredName.split(",");
						Disjunction disj = Restrictions.disjunction();
						for (String val : vals)
							disj.add(Restrictions.eq("mp2.id", Integer.valueOf(val.substring(1))));
						criteria.add(disj);
					}
				}
				else if (enteredName.matches("M[0-9]*\\*")) // MXXXXX* (mapping2 with stereoisomers)
				{
					Mapping2 mp2 = (Mapping2) Globals.session().get(Mapping2.class, Integer.valueOf(enteredName.substring(1, enteredName.length() - 1)));
					criteria.add(Restrictions.eq("mp1.id", mp2.mapping1.id));
				}
				else if (enteredName.matches("R[0-9]*")) // RXXXX - record ID
					criteria.add(Restrictions.eq("id", Long.valueOf(enteredName.substring(1))));
				else if (enteredName.matches("R[0-9]*\\*")) // RXXXX* - find similar records
				{
					ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, Long.valueOf(enteredName.substring(1, enteredName.length() - 1)));
					criteria.add(Restrictions.eq("mp1.id", ep.molecule.mapping1.id));
					criteria.add(Restrictions.eq("property", ep.property));
				}
				else if (enteredName.matches("Q[0-9]*")) // internal depiction id
					criteria.add(Restrictions.eq("id", Long.valueOf(enteredName.substring(1))));
				else if (enteredName.matches("[0-9]*"))
					// Regardless of stereochemistry
					criteria.add(Restrictions.eq("mp1.id", filters.getLong("name")));
				else if (enteredName.matches("ID2=[0-9]*"))
					// Considering stereochemistry
					criteria.add(Restrictions.eq("mp2.id", Integer.valueOf(enteredName.substring(4))));
				else if (enteredName.matches("ID1=[0-9]*"))
					// Considering stereochemistry
					criteria.add(Restrictions.eq("mp1.id", Integer.valueOf(enteredName.substring(4))));
				else if (enteredName.matches("[A-Z]{14}"))
				{
					// Inchie-1
					criteria.add(Restrictions.eq("mp1.inchi1", enteredName));
				}
				else if (enteredName.matches("[A-Z]{14}-[A-Z]{10}"))
					criteria.add(Restrictions.eq("mp2.inchi2", enteredName));
				else
				{
					// Other than that - name
					distinctQueryRequired = true;
					String checkIfCasrn = CASRN.checkCasrnSyntax(enteredName);
					if (null != checkIfCasrn)
						enteredName = checkIfCasrn;

					criteria
					.createCriteria("moleculenames")
					//						.createCriteria("name")
					.add(Restrictions.like("name", enteredName+"%"));
				}
			}

			if (filters.has("approval-status"))
			{
				String value = filters.get("approval-status");
				if ("only-approved".equals(value) && filters.has("introducer"))
					criteria.add(Restrictions.eq("approved", true));
				if ("only-awaiting-approval".equals(value))
				{
					criteria.add(Restrictions.eq("approved", false));
					criteria.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE));
					criteria.add(Restrictions.eq("rejected", false));
				}

				if ("rejected".equals(value)) {
					criteria.add(Restrictions.eq("approved", false));
					criteria.add(Restrictions.eq("rejected", true));
					criteria.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE));	
				}
			}

			if (filters.has("visibility"))
			{
				String value = filters.get("visibility");
				if ("public".equals(value))
					criteria.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE));
				else if ("private".equals(value))
					criteria.add(Restrictions.eq("rights", Globals.RIGHTS_NONE));
			}

			if (filters.has("basket-select") && !filters.get("basket-select").startsWith("-"))
			{
				// Normal basket (union)
				Disjunction junction = Restrictions.disjunction();
				if (filters.get("basket-select").equals("any"))
					distinctQueryRequired = true;
				else
				{
					String ids[] = filters.get("basket-select").split(",");
					if (ids.length > 1)
						distinctQueryRequired = true;
					for (String id : ids)
					{
						Basket basket = null;
						try
						{
							long basketId = Long.valueOf(id);
							basket = Basket.getBasket(Globals.userSession(), basketId);
							if(basket == null)throw new Exception("not found " +basketId );
						} catch (Exception e)
						{
							String basketName = id;
							basket = Basket.getBasket(Globals.userSession(), basketName);
						}
						junction.add(Restrictions.eq("b.id", basket == null ? -1: basket.id));
					}
				}
				criteria.add(junction);

				criteria.createAlias("basketEntries", "be");
				criteria.createAlias("be.basket", "b");

				// Show only the excluded records
				if (assertParam("basket-excluded"))
					criteria.add(Restrictions.eq("be.exclude", Boolean.TRUE));

				// And here we show the ercords excluded for a particular model
				if (assertParam(QSPRConstants.EXCLUDED))
				{
					Disjunction excluded = Restrictions.disjunction();

					Long id = getLongParam(QSPRConstants.EXCLUDED); //Model ID for which we're looking for excluded. If negative, means EVERYTHING but excluded
					Model m = (Model)Globals.session().get(Model.class, Math.abs(id));

					if (m.microattachment.getObject().excludedBasketEntries.size() > 0)
						excluded.add(Restrictions.in("be.id", m.microattachment.getObject().excludedBasketEntries));
					else
						excluded.add(Restrictions.sqlRestriction("1 = 0"));

					excluded.add(Restrictions.eq("be.exclude", true));

					if (id > 0)
						criteria.add(excluded);
					else
						criteria.add(Restrictions.not(excluded));
				}
			}

			if (filters.has("modelerrors"))
			{
				// Show records that are errors in given model
				// Just use by-list filtering, similar to what we do for selected items
				// Midnighter

				Model model = (Model) Globals.session().get(Model.class, getLongParam("modelerrors"));
				if (model == null)
					model = (Model) Globals.getSessionAttribute(SessionVariable.MODEL);

				List<Long> errorIds = new ArrayList<Long>();
				for (ModelMapping mm : model.modelMappings) 
					errorIds.addAll(ModelDotService.getErrorRecordsIds(mm, assertParam("recalculated"), getParam("modelError")));

				if (errorIds == null || errorIds.size() == 0)
					criteria.add(Restrictions.sqlRestriction("1 = 0"));
				else
					criteria.add(Restrictions.in("id", errorIds));
			}

			if (filters.has("duplicates"))
				if (filters.has("article") || (filters.has("basket-select") && !filters.get("basket-select").startsWith("-")))
				{
					// Draft version of dublicates search. Implemented using subquery
					// Dublicates by property+molecule are considered.
					// If article or basket filters are activated, we restrict dublicate search to appropriate context
					// Midnighter

					//criteria.createAlias("mol.mapping2", "mp2");

					DetachedCriteria dublicatesCount = 
							DetachedCriteria.forClass(ExperimentalProperty.class, "d");
					dublicatesCount.createAlias("molecule", "dmol");
					dublicatesCount.createAlias("property", "dp");

					// Duplicates: Same molecule
					if (filters.has("nostereo"))
					{
						dublicatesCount.createAlias("dmol.mapping1", "dmp1");
						dublicatesCount.add(Restrictions.eqProperty("mp1.id", "dmp1.id"));					
					} 
					else
					{
						dublicatesCount.createAlias("dmol.mapping2", "dmp2");
						dublicatesCount.add(Restrictions.eqProperty("mp2.id", "dmp2.id"));
					}

					dublicatesCount
					.add(Restrictions.eqProperty("p.id", "dp.id"))
					.add(Restrictions.isNull("deleted"))
					.setProjection(Projections.countDistinct("d.id"));

					if (filters.has("article"))
						dublicatesCount.add(Restrictions.eqProperty("this.article", "d.article"));
					if (filters.has("basket-select") && !filters.get("basket-select").startsWith("-"))
					{
						// A normal basket
						Disjunction junction = Restrictions.disjunction();
						String ids[] = filters.get("basket-select").split(",");
						for (String id : ids) {
							Basket basket = Basket.getBasket(Globals.userSession(), Long.valueOf(id));
							junction.add(Restrictions.eq("id", basket.id));
						}
						dublicatesCount
						.createCriteria("basketEntries")
						.createCriteria("basket")
						.add(junction);
					}

					criteria.add(Subqueries.leAll(2L, dublicatesCount));
				}

			if (filters.has("compare-baskets") && filters.get("basket-select").contains(","))
			{
				DetachedCriteria duplicatesCount = DetachedCriteria.forClass(BasketEntry.class, "be");
				duplicatesCount.createAlias("be.ep", "dep");
				duplicatesCount.createAlias("dep.molecule", "dmol");
				duplicatesCount.createAlias("dmol.mapping2", "dmp2");
				duplicatesCount.createAlias("dep.property", "dp");

				duplicatesCount
				.add(Restrictions.eqProperty("mp2.id", "dmp2.id"))
				.add(Restrictions.eqProperty("p.id", "dp.id"))
				.setProjection(Projections.countDistinct("basket.id"));

				Disjunction junction = Restrictions.disjunction();
				String ids[] = filters.get("basket-select").split(",");
				for (String id : ids) {
					Basket basket = Basket.getBasket(Globals.userSession(), Long.valueOf(id));
					logger.info(basket);
					junction.add(Restrictions.eq("be.basket", basket));
				}

				duplicatesCount.add(junction);

				criteria.add(Subqueries.leAll(Long.valueOf(ids.length), duplicatesCount));
			}

			if (ssFilter.structure != null)
				criteria.setMaxResults(1000);

			Order order;	

			if (!globalSSSearch)
				if (filters.has("sortby"))
				{
					if (filters.get("sortby").equals("molprop"))
					{
						order = (filters.has("order")) ? Order.asc("mp1.id") : Order.desc("mp1.id");
						wrapper.order.add(order);
						order = (filters.has("order")) ? Order.asc("mp2.id") : Order.desc("mp2.id");
						wrapper.order.add(order);					
						order = (filters.has("order")) ? Order.asc("property") : Order.desc("property");
						wrapper.order.add(order);
					} else
					{
						order = (filters.has("order")) ? Order.asc(filters.get("sortby")) : Order.desc(filters.get("sortby"));
						wrapper.order.add(order);
						// put id as second order criteria if records are the same by the first
						order = (filters.has("order")) ? Order.asc("id") : Order.desc("id");
						wrapper.order.add(order);
					}
				} else
				{
					order = (filters.has("order")) ? Order.asc("id") : Order.desc("id");
					wrapper.order.add(order);				
				}
		}

		if (filters.has("trash"))
		{
			criteria.add(Restrictions.isNotNull("deleted"));
			if (user != null)
			{
				if (!filters.has("globaltrash") || !Globals.userSession().user.isSuperUser())
				{
					criteria.createAlias("introducer", "int", JoinType.LEFT_OUTER_JOIN);
					criteria.add(Restrictions.eq("int.id", user.id));
				}
			}
			else
				criteria.add(Restrictions.isNull("introducer"));
		}
		else
		{
			criteria.add(Restrictions.isNull("deleted")); 
			boolean showAwaitingApproval = assertParam("awaiting-approval");
			if (showAwaitingApproval)
			{
				criteria.add(Restrictions.eq("approved", false));
				criteria.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE));
				wrapper.createAlias("property", "p");
				if (!assertParam("all-moderators"))
				{
					Disjunction disjunction = Restrictions.disjunction();
					disjunction.add(Restrictions.eq("p.moderator", Globals.userSession().user)); // either moderated by me
					if (assertParam("unmoderated"))
						disjunction.add(Restrictions.isNull("p.moderator"));// or with no moderators yet
					criteria.add(disjunction);
				}
			}

			boolean showUnappproved = !"only-approved".equals(filters.get("approval-status"));
			ExperimentalProperty.addAccessRestrictions(criteria, rights, null, wrapper.aliases.containsKey("b"), showUnappproved);
		}

		distinctRequired.value = distinctQueryRequired;

		return wrapper;
	}

	protected void formFilter(String key, String value, WebFilters filters)
	{
		if (key.equals("basket-select"))
		{
			Basket basket = null;
			if ("selected".equals(value))
				basket = Globals.selectionBasket(false);
			else if (!value.contains(",") && !value.equals("any"))
			{	
				try
				{
					basket = Basket.getBasket(Globals.userSession(), getLongParam("basket-select"));
				} catch (NumberFormatException e)
				{
					basket = Basket.getBasket(Globals.userSession(), getParam("basket-select"));
				}
			}

			if (basket != null)
				filters.addFilter(key, "" + basket.id, basket.name + " (" + basket.session.user + ")");
			else
				filters.addFilter(key, value, "");
		}
		else if (key.equals("basket"))
		{
			Criteria criteria = Globals.session().createCriteria(Basket.class);
			criteria.add(Restrictions.eq("name", value));
			Session session = Globals.userSession();
			if (session.user != null)
				criteria.add(Restrictions.eq("user", session.user));
			else
				criteria.add(Restrictions.eq("session", session));
			List<Basket> baskets = criteria.list();	
			if (baskets.size() != 0)
			{
				Basket basket = baskets.get(0);
				filters.addFilterOverride("basket-select", basket.id.toString(), "Basket filter");
			}
		} else

			if (key.equals("article"))
			{
				Article article = (Article) Globals.session().get(Article.class, Long.valueOf(value));
				filters.addFilter("article", (Long.valueOf(article.id)).toString(), article.getTitle());
			} else
				if (key.equals("property"))
				{
					if (!value.contains(","))
					{
						Property property = (Property) Globals.session().get(Property.class, Long.valueOf(value));
						filters.addFilter("property", (Long.valueOf(property.id)).toString(), property.getName());
					}
				} else
					if (key.equals("similarto"))
					{
						ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, Long.valueOf(value));
						if (ep == null)
							throw new UserFriendlyException("Property not found");
						filters.addFilter("property", (Long.valueOf(ep.property.id)).toString(), ep.property.getName());

						if (assertParam("stereochemistry"))
							filters.addFilter("name", "M" + ep.molecule.mapping2.id.toString(), "");
						else
							filters.addFilter("name", "M" + ep.molecule.mapping2.id.toString() + "*", "");


						//filters.addFilter("molid", ep.molecule.id.toString(), "");
						if (ep.property.isNumeric())
							filters.addFilter("canonicalValue", ep.canonicalValue.toString(), "");
					} else
						if (key.equals("similarmol"))
						{
							if (Long.parseLong(value) > 0)
							{
								Molecule molecule = Repository.molecule.getMolecule(Long.valueOf(value));
								if (assertParam("stereochemistry"))
									filters.addFilter("name", "ID2=" + molecule.mapping2.id.toString(), "");
								else
									filters.addFilter("name", molecule.mapping1.id.toString(), "");
							} else {
								// TODO: think about similar depiction for pubchem mols
								filters.addFilter("name", "-999", "");
							}
						} else
							super.formFilter(key, value, filters);
	}

	public ModelAndView list(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		WebList list = new WebList();
		WebFilters filters = formFilters(req);

		RequestStatistics stats = new RequestStatistics();
		stats.start();
		Mutable<Boolean> distinctRequired = new Mutable<Boolean>();

		setStatus("Preparing the query...");
		CriteriaWrapper criteriaWrapper = createMainCriteria(filters, Globals.RIGHTS_FREELY_AVAILABLE, Globals.session().createCriteria(ExperimentalProperty.class), distinctRequired);
		Criteria criteria = criteriaWrapper.criteria;
		list.setCriteriaWrapper(criteriaWrapper);

		if (filters.has("basket-select") && !filters.get("basket-select").equals("-1")){
			list.maxAllowedItemsForSorting = 1000000; // If we select a basket, allow sorting for even very large result sets
			Globals.setSessionAttribute(SessionVariable.BASKET_SELECT, filters.getLong("basket-select"));
		}
		else
			Globals.setSessionAttribute(SessionVariable.BASKET_SELECT, null);

		setStatus("Quering results...");
		logger.info("Criteria created: " + stats.current());

		if (distinctRequired.value)
			list
			.useEntity(ExperimentalProperty.class)
			.loadDistinctFromCriteria(criteria, getPageNum(), getPageSize(5));
		else
			list
			.useEntity(ExperimentalProperty.class)
			.loadFromCriteria(criteria, getPageNum(), getPageSize(5));

		stats.stop();
		logger.info("List: "+stats);

		if (assertParam("rights"))
			if (Globals.userSession().user != null)
				Globals.userSession().user.defaultRights = Integer.valueOf(req.getParameter("rights"));

		BrowserModel browserModel = new BrowserModel().setFilters(filters);
		if (Globals.userSession().user != null)
			browserModel.selectionSize = (Integer)Globals.session().createSQLQuery("select count(*) c from BasketEntry natural left join ExperimentalProperty ep left join Basket b using (basket_id) where b.user_id="+Globals.userSession().user.id+" and b.basket_type=1 and ep.deleted is null")
			.addScalar("c", IntegerType.INSTANCE)
			.uniqueResult();

		return browserModel.setObject(list).getModelAndView();
	}

	public ModelAndView getduplicate(HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		Long id = getLongParam("id");
		ExperimentalProperty ep = (ExperimentalProperty)Globals.session().get(ExperimentalProperty.class, id);

		Pattern regex = Pattern.compile(".*R([0-9]+).*",Pattern.DOTALL);

		if (ep.errorComment == null || !ep.errorComment.startsWith("Duplicate"))
			return new WebModel(null).getModelAndView();

		Matcher regexMatcher = regex.matcher(ep.errorComment);

		if (!regexMatcher.matches())
			return new WebModel(null).getModelAndView();

		id = Long.valueOf(regexMatcher.group(1));
		ep = (ExperimentalProperty)Globals.session().get(ExperimentalProperty.class, id);

		return new WebModel(ep).getModelAndView();
	}

	private void executeDelete(ExperimentalProperty ep, HttpServletRequest request) throws Exception
	{
		boolean permanentDeletion = assertParam("trash");
		if (permanentDeletion)
		{
			// Only introducer or superuser can delete record permanently
			if (!ep.isDeleted())
				return;
			if (ep.introducer != null && !ep.introducer.equals(Globals.userSession().user) && !Globals.userSession().user.isSuperUser())
				return;
		}
		else
			AccessChecker.requestModificationPermission(ep);

		// 2. Delete all my basket entries connected with this record
		List<BasketEntry> bes = Globals.session()
				.createCriteria(BasketEntry.class)
				.add(Restrictions.eq("ep", ep))
				.list();

		for (BasketEntry basketEntry : bes)
		{
			basketEntry.basket.cachedCount = null;
			Globals.session().delete(basketEntry);
		}

		// 3. Unbind records referencing this one
		ep.colorednames = null;
		List<ExperimentalProperty> referencingEps = (List<ExperimentalProperty>) Globals.session().createCriteria(ExperimentalProperty.class)
				.add(Restrictions.eq("connectedProperty", ep)).list();
		for (ExperimentalProperty experimentalProperty : referencingEps) {
			experimentalProperty.connectedProperty = null;
		}
		Globals.session().flush();

		// 4. Finally delete the record itself
		if (permanentDeletion)
			Globals.session().delete(ep);
		else
		{
			// Move the record to the trash of the introducer
			ep.deleted = new Timestamp(Calendar.getInstance().getTimeInMillis());
			ep.owner = Globals.userSession().user;
			if (ep.introducer == null)
				ep.introducer = Globals.userSession().user;

			ep.updateHash();
			Globals.session().saveOrUpdate(ep);
		}

		// 5. and also remove it from selectionList
		Globals.selectionBasket(assertParam("trash")).removeEntry(ep.id);
		//Globals.userSession().selectionList.remove(ep.id);

	}

	private void executeAction(ExperimentalProperty ep, HttpServletRequest request) throws Exception
	{
		String action = request.getParameter("action");
		if (action.equals("delete") || action.equals("deleteselected"))
			executeDelete(ep, request);
		else if (action.equals("removeselect") || action.equals("addselect") || action.equals("toggleselect") || action.equals("selectpage") || action.equals("cleanbasket"))
		{
			//Set<Long> selectionList = Globals.userSession().selectionList;

			if (action.equals("removeselect"))
			{
				Globals.selectionBasket(assertParam("trash")).removeEntry(ep.id);
				//				selectionList.remove(ep.id);
			}
			else if (action.equals("addselect") || action.equals("selectpage"))
			{
				Globals.selectionBasket(assertParam("trash")).addEntry(ep.id);
				//				selectionList.add(ep.id);
			}
			else if(action.equals("cleanbasket"))
			{
				String hash = ""+ep.property.id;
				if(request.getParameter("withStreo").equals("true"))
					hash += "_"+ep.molecule.mapping2.id;
				else
					hash += "_"+ep.molecule.mapping1.id;

				if(request.getParameter("withValue").equals("true"))
					hash += "_"+ep.value;

				//create hash
				String md5 = OCHEMUtils.getMD5(hash);
				Long dubRowNum = md5Hashes.get(md5);

				if(dubRowNum != null)
				{
					if(ep.connectedProperty != null && ep.id == ep.connectedProperty.id)
					{
						Globals.selectionBasket(assertParam("trash")).addEntry(dubRowNum);
						//						selectionList.add(dubRowNum);
					}
					else
					{
						Globals.selectionBasket(assertParam("trash")).addEntry(ep.id);
						//						selectionList.add(ep.id);
					}
				}
				else
					md5Hashes.put(md5, ep.id);

			}
			else
			{
				// Toggle selection
				//				if (selectionList.contains(ep.id))
				//					selectionList.remove(ep.id);
				//				else
				//					selectionList.add(ep.id);
				if (Globals.selectionBasket(assertParam("trash")).containsEntry(ep.id))
					Globals.selectionBasket(assertParam("trash")).removeEntry(ep.id);
				else
					Globals.selectionBasket(assertParam("trash")).addEntry(ep.id);
			}

		}
		else if ((action.equals("removebasket"))  || (action.equals("addbasket")))
		{
			Session currentSession = Globals.userSession();
			Basket currentBasket = Basket.getBasket(currentSession, request.getParameter("basket-name"), false);
			currentBasket.cachedCount = null;

			Criteria c = Globals.session()
					.createCriteria(BasketEntry.class)
					.add(Restrictions.eq("basket", currentBasket));
			if(assertParam("new-mol"))
			{
				c.createAlias("ep.molecule", "mol")
				.add(Restrictions.eq("mol.id", ep.molecule.id));
			}
			else if(assertParam("ignore-stereo"))
			{
				c.createAlias("ep.molecule", "molecule")
				.createAlias("molecule.mapping1", "mp1")
				.add(Restrictions.eq("mp1.id", ep.molecule.mapping1.id));
			}
			else
				c.add(Restrictions.eq("ep", ep));

			List<BasketEntry> entries = c.list();

			if (entries.size() > 0)
			{
				if (!action.equals("addbasket"))
					for (BasketEntry entry : entries) 
						Globals.session().delete(entry);
			}
			else 
				if (!action.equals("removebasket"))
				{
					BasketEntry be = new BasketEntry();
					be.ep = ep;
					be.basket = currentBasket;
					currentSession.mod_time = new Timestamp(Calendar.getInstance().getTimeInMillis());
					Globals.session().save(currentSession);
					Globals.session().save(be);
				}
		}
		else if (action.equals("restore") || action.equals("restoreselected"))
		{
			// Restore record from trash
			if (ep.introducer != Globals.userSession().user && !Globals.userSession().user.isSuperUser())
				throw new UserFriendlyException("Not sufficient privileges to access model of user: "
						+ Globals.userSession().user.login);
			ep.deleted = null;
			ep.updateHash();
			//add handler here
			if (ep.hasConflicts())
				if (ep.isPublic() && !ep.duplicate.isPublic())
				{
					// Public invalidates hidden
					ep.duplicate.ep_status = 0;
					ep.duplicate.updateHash();
					Globals.session().saveOrUpdate(ep.duplicate);
					Globals.session().flush();
				}
				else 
					throw new DublicateRestoreException(ep.duplicate);
			Globals.session().saveOrUpdate(ep);
		}
		else if (action.equals("addtag"))
		{
			// Add tag to compound
			Tag tag = (Tag) Globals.session().get(Tag.class, getLongParam("tagid"));
			if (!ep.molecule.isEmptyMolecule() && !ep.molecule.mapping1.tags.contains(tag))
			{
				ep.molecule.mapping1.tags.add(tag);
				Globals.session().saveOrUpdate(ep.molecule.mapping1);
			}
		}
		else if (action.equals("removetag"))
		{
			// Remove tag from compound
			Tag tag = (Tag) Globals.session().get(Tag.class, getLongParam("tagid"));
			if (!ep.molecule.isEmptyMolecule() && ep.molecule.mapping1.tags.contains(tag))
			{
				ep.molecule.mapping1.tags.remove(tag);
				Globals.session().saveOrUpdate(ep.molecule.mapping1);
			}
		}
		else if (action.equals("publish"))
		{
			AccessChecker.requestModificationPermission(ep);
			ep.rights = Globals.RIGHTS_FREELY_AVAILABLE;
			ep.rejected = false;
			ep.updateHash();
			if (!ep.hasConflicts())
			{
				Globals.session().saveOrUpdate(ep);
			}
			else
				throw new UserFriendlyException("You cant publish this record, because it will create dublicated entry in public domain");
		}
		else if (action.equals("approve") || action.equals("disapprove") || action.equals("unapprove"))
		{
			//			ep.approve();
			Session currentSession = Globals.userSession();
			Basket currentBasket = Basket.getBasket(currentSession, "Basket of reviewed records", false);
			BasketEntry be = new BasketEntry();
			be.ep = ep;
			be.basket = currentBasket;
			currentSession.mod_time = new Timestamp(Calendar.getInstance().getTimeInMillis());
			Globals.session().save(currentSession);
			Globals.session().save(be);
		}
		//		else if ()
		//		{
		//			// Basically, "unpublish" the record
		//			if (ep.property.moderator != null && !ep.property.moderator.equals(Globals.userSession().user))
		//				throw new UserFriendlyException("You are not allowed to approve a record for the property " + ep.property.getName() + ", which is moderated by " + ep.property.moderator.login);
		//			ep.approved = false;
		//			ep.rejected = true;
		//			//ep.rights = Globals.RIGHTS_NONE;
		//		}
		else 
			throw new Exception("Unknown action - "+action);
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		String action = request.getParameter("action");
		if (action.equals("clearselection"))
		{
			Globals.selectionBasket(false).removeEntries();
			return new WebModel().getModelAndView();
		}

		WebFilters filters = formFilters(request);

		if ("approve".equals(action) || "disapprove".equals(action) || "unapprove".equals(action))
		{
			new DataApprovalEmailNotifier();
			Session currentSession = Globals.userSession();
			Basket currentBasket = Basket.getBasket(currentSession, "Basket of reviewed records", false);
			currentBasket.cachedCount = null;
			currentBasket.removeEntries();
			Globals.session().flush();
			Globals.session().saveOrUpdate(currentBasket);
		}

		if("similar-record".equals(action))
		{

			ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, getLongParam("id"));
			Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class)
					.createAlias("molecule", "mol")
					.add(Restrictions.eq("article", ep.article))
					.add(Restrictions.ne("id", ep.id))
					.add(Restrictions.isNull("deleted"));
			if(ep.moleculenames.size() > 0)
			{
				criteria.createCriteria("moleculenames")
				.add(Restrictions.like("name", "%"+ep.moleculenames.get(0)+"%"));
			}
			criteria
			.add(Restrictions.eq("mol.mapping1", ep.molecule.mapping1));
			WebList list = new WebList();
			list.useEntity(ExperimentalProperty.class).loadDistinctFromCriteria(criteria, getPageNum(), getPageSize(5));
			return new WebModel(list).setTemplate("similar-record-browser").getModelAndView();
		}

		int maxAllowed = Globals.userSession().getLimit();
		if ("removebasket".equals(action))
		{	
			filters.addFilter("apply-basket", "1", "");
			//maxAllowed = 10000;
		}

		//remove duplicates from basket currently duplicate criteria is based on mapping_id1, property and value
		if("cleanbasket".equals(action))
			md5Hashes.clear();

		if (filters.in("action", new String[]{"deleteselected", "delete"}) && !assertParam("trash"))
		{
			// Get introducers of the data to delete
			List<Object[]> affectedUsersRows = null; 
			{
				Mutable<Boolean> distinctRequired = new Mutable<Boolean>();
				Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class);
				createMainCriteria(filters, Globals.RIGHTS_FREELY_AVAILABLE, criteria, distinctRequired);
				criteria.setProjection(Projections.projectionList().add(Projections.groupProperty("introducer")).add(Projections.count("id")));
				affectedUsersRows = criteria.list();
			}

			// Get baskets affected by deletion (wolfram->)			
			// e.g. SELECT distinct basket_id FROM qspr.ExperimentalProperty E join BasketEntry F on (E.exp_property_id = F.exp_property_id) where E.exp_property_id < 10
			// or SELECT distinct basket_id, introducer_id FROM qspr.ExperimentalProperty E join BasketEntry F on (E.exp_property_id = F.exp_property_id) where E.exp_property_id<30 order by introducer_id						
			List<Basket> list2 = null; 
			{
				// Case one : if we have already joined with basketentrys, two queries must be issued. If not, we can add the join() to the main criterion and execute only one query ("else")
				if (filters.has("basket-select") && !filters.get("basket-select").startsWith("-"))
				{
					// (Two step process necesssary) - two joins with basketentry table :(
					// First join()
					List<ExperimentalProperty> lep;
					{
						Mutable<Boolean> distinctRequired = new Mutable<Boolean>();
						Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class);
						createMainCriteria(filters, Globals.RIGHTS_FREELY_AVAILABLE, criteria, distinctRequired);
						criteria.setProjection(Projections.id());
						lep = criteria.list();						
					}
					// Second join()
					Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class);
					criteria.add(Restrictions.in("id", lep.toArray()));

					criteria.createAlias("basketEntries", "be");
					criteria.setProjection(Projections.projectionList().add(Projections.groupProperty("be.basket")));
					list2 = criteria.list();

				} else {
					Mutable<Boolean> distinctRequired = new Mutable<Boolean>();
					Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class);
					createMainCriteria(filters, Globals.RIGHTS_FREELY_AVAILABLE, criteria, distinctRequired);
					criteria.createAlias("basketEntries", "be");					
					criteria.setProjection(Projections.projectionList().add(Projections.groupProperty("be.basket")));
					list2 = criteria.list();
				}
			}			

			// We have the baskets affected in list2 - BasketEntries contain Name of the Basket and User Id
			// where basket_type=0 (1 means 'selected records')

			// find affected users that get mail message if del() is executed:
			List<User> affectedBasketsOwners = new ArrayList<User>();
			StringBuilder basketDeletionSummaryMessage = new StringBuilder();
			{
				for (Basket basket : list2)
					if (basket.user != null) // We care about registered users only
					{
						if (basket.basketType == 0L) 
						{ // User basket
							basketDeletionSummaryMessage.append("Basket \""+basket.name+"\" from user \""+basket.user.login+"\"\n");
							if (!affectedBasketsOwners.contains(basket.user)) 
								affectedBasketsOwners.add(basket.user); // This user is affected.
						}					
					}
			}			

			boolean messageRequired = false;
			List<User> introducers = new ArrayList<User>();
			int numOfRecords = 0;
			String deletionSummaryMessage = "";
			for (Object[] objects : affectedUsersRows)
			{
				if (objects[0] != null)
					if (!((User)objects[0]).equals(Globals.userSession().user))
					{
						messageRequired = true;
						introducers.add((User)objects[0]);
						deletionSummaryMessage += objects[0]+": "+objects[1]+"\n";
					}
					else
						deletionSummaryMessage += "your own records: "+objects[1]+"\n";
				else
					deletionSummaryMessage += "anonymous records: "+objects[1]+"\n";
				numOfRecords += (Long)objects[1];
			}

			if (numOfRecords == 0)
				throw new UserFriendlyException("You did not select any records");

			if (!assertParam("deletion-confirmation"))
			{
				// Send selected data overview back to the user, and ask for his confirmation
				String msg = "You are about to *PERMANENTLY* delete "+numOfRecords+" records from the database. Are you sure that you want to do it?";
				if (basketDeletionSummaryMessage.length()>0) {
					msg += "\n\nThis delete will affect the following user baskets : \n"+basketDeletionSummaryMessage+"The records will be removed from these baskets.\n\n";	
					//msg += "\nSome of these records are currently contained in user baskets : \n"+basketDeletionSummaryMessage+"\nThe records will be removed from these baskets.";					
				}
				if (messageRequired)
					return new WebModel(new Alert(msg + "Some of them belong to other user(s):\n\n"
							+deletionSummaryMessage+"\nPlease provide a reason for deletion for these users", "provide-message")).getModelAndView();
				else
					return new WebModel(new Alert(msg + "Are you sure you really want to delete all these records?")).getModelAndView();
			}
			else 
			{
				// So we have confirmation. Send messages to introducers of data to be deleted
				for (Object[] objects : affectedUsersRows) 
				{
					if (objects[0] != null)
						if (!((User)objects[0]).equals(Globals.userSession().user))
						{
							User _user = Globals.userSession().user;
							if(_user == null)
								throw new UserFriendlyException("Not permitted, delete has been ignored");
							// Messages to data introducers
							String subject = ""+objects[1]+" of your records have been deleted by "+Globals.userSession().user+ " at "+Globals.userSession().time.toString()+"\n\n";
							String body = subject+".\n You can find and review them in your trash.\n\n\n"+"\n\n"+"Reason for deletion:\n" + getParam("deletion-confirmation")+"\n\n";
							dialogueService.sendMessage(_user, (User)objects[0], subject, body, true);								
						}
				} 

				// Send messages to users that have records in baskets
				for (User user : affectedBasketsOwners) 
				{
					if (user.equals(Globals.userSession().user))
						continue;

					String affectedBasketsMsg = "";
					for (Basket basket : list2) 
						if (basket.user == user)  // This basket belongs to the current user.
							affectedBasketsMsg += ","+basket.name;
					if (affectedBasketsMsg.length() > 0) 
						affectedBasketsMsg = affectedBasketsMsg.substring(1);

					String subject = "Records that are used in your baskets have been deleted by "+Globals.userSession().user+ " at "+Globals.userSession().time.toString()+".\n\n";
					String body = subject + "\nAffected baskets: "+"\n\n"+affectedBasketsMsg+"."+"\n\n";
					body += "\n\nReason for deletion:\n"+"\n\n" + getParam("deletion-confirmation")+"."+"\n\n";
					dialogueService.sendMessage(Globals.userSession().user, user, subject, body, true);								
				} 

			}
		}

		Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class);
		Mutable<Boolean> distinctRequired = new Mutable<Boolean>();
		createMainCriteria(filters, Globals.RIGHTS_FREELY_AVAILABLE, criteria, distinctRequired);

		ExperimentalProperty returnedObject = null;

		if (filters.has("id") && filters.get("id").equals("-1"))
		{
			returnedObject = new ExperimentalProperty();

			if (filters.has("hide") && filters.get("hide").equals("true"))
				returnedObject.rights = Globals.RIGHTS_NONE;

			executeAction((ExperimentalProperty) returnedObject, request);
		}
		else
		{
			WebList l = new WebList();

			// 1. Select all IDs of matched records
			// 2. Fetch objects for these IDs, by 50 records. 
			//    Restart transaction after each 50-records job

			l.useEntity(ExperimentalProperty.class);
			if (distinctRequired.value)
				l.useProjection(Projections.distinct(Projections.property("id")));
			else
				l.useProjection(Projections.property("id"));
			logger.info("Starting the query...");
			l.loadFromCriteria(criteria, 1, maxAllowed);

			List<Object> list = l.list;
			if (l.size > maxAllowed)
				throw new UserFriendlyException("You tried to work with " + l.size + " records. The maximum allowed number is currently "+maxAllowed+" records");

			// Do not add a records to a basket twice
			List<Long> currentBasketRecords = null;

			if ("addbasket".equals(filters.get("action")))
			{
				Basket currentBasket = Basket.getBasket(Globals.userSession(), request.getParameter("basket-name"), false);
				if (filters.has("clearbasket"))
				{
					currentBasket.removeEntries();
					Globals.session().saveOrUpdate(currentBasket);
					Globals.session().flush();
				}
				currentBasketRecords = Globals.session().createCriteria(BasketEntry.class).createAlias("ep", "ep_a").add(Restrictions.eq("basket", currentBasket)).setProjection(Projections.groupProperty("ep_a.id")).list();
				list.removeAll(currentBasketRecords);
			}

			Iterator<Object> ieps = list.iterator();
			int counter = 0;
			int batchSize = isBatchAction(action) ? 100000 : 50;
			// A loop processing records in batches
			while(ieps.hasNext())
			{
				int i = 0;
				List<Long> ids = new ArrayList<Long>();
				while (ieps.hasNext() && i++ < batchSize)
					ids.add((Long)ieps.next());

				counter += ids.size();

				Criteria criteria1 = Globals.session().createCriteria(ExperimentalProperty.class)
						.add(Restrictions.in("id", ids));

				if("cleanbasket".equals(filters.get("action")))
				{
					criteria1.createAlias("molecule", "mol")
					.createAlias("mol.mapping1", "mp1")
					.addOrder(Order.asc("mp1.id"))
					.createAlias("article", "art")
					.addOrder(Order.asc("art.publicationDate"));
				}

				if (isBatchAction(action))
				{
					logger.info("Executing " + filters.get("action") + " for " + ids.size() + " records");
					if ("addbasket".equals(filters.get("action")) || "removebasket".equals(filters.get("action")))
					{
						// For sake of efficienty, here goes a direct batch SQL query
						Basket currentBasket = Basket.getBasket(Globals.userSession(), request.getParameter("basket-name"), false);
						if ("addbasket".equals(filters.get("action")))
						{
							int mol_option = (filters.has("mol-option")) ? filters.getInteger("mol-option") : 0;
							currentBasket.addEntries(ids, mol_option);
						} else
							currentBasket.removeEntries(ids);
					}

					if ("addselect".equals(filters.get("action")))
					{
						Globals.selectionBasket(assertParam("trash")).addEntries(ids);
						//Globals.userSession().selectionList.addAll(ids);
					}

					if ("removeselect".equals(filters.get("action")))
					{
						Globals.selectionBasket(assertParam("trash")).removeEntries(ids);
						//Globals.userSession().selectionList.removeAll(ids);
					}

					if ("publish".equals(filters.get("action")) || ACTION_PUBLISH_FREE.equals(filters.get("action")))
					{
						if (Globals.userSession().user == null)
							throw new UserFriendlyException("Guest users cannot publish data");

						//List<Long> l= Globals.session().createQuery("select publicId from Model where id=:ids").setParameter("ids",id)

						ExperimentalProperty.checkIfSetIsPublishable(ids);

						// Modify access rights
						Globals.session().createQuery("update ExperimentalProperty set rights=" + Globals.RIGHTS_FREELY_AVAILABLE + 
								" where id in (:ids) and introducer=:introducer and rights < " + Globals.RIGHTS_FREELY_AVAILABLE)
						.setParameterList("ids", ids)
						.setParameter("introducer", Globals.userSession().user)
						.executeUpdate();
					}
				}
				else
				{
					// Execute the action records by record
					List<ExperimentalProperty> eps = criteria1.list();
					for (ExperimentalProperty ep : eps) {
						returnedObject = ep;
						logger.info("[EBC] Executing "+request.getParameter("action")+" on record "+ep.id);
						executeAction((ExperimentalProperty)ep, request);
					}
				}

				logger.info(getParam("action") + ": " + counter + " out of "+list.size()+" records processed");

				Globals.restartAllTransactions(true);
			}

			if ("addbasket".equals(filters.get("action")))
				EventFactory.document("Add to basket", new BasketAction(Basket.getBasket(Globals.userSession(), request.getParameter("basket-name")), list.size(), BasketActionType.ADDENTRY));
			if ("removebasket".equals(filters.get("action")))
				EventFactory.document("Add to basket", new BasketAction(Basket.getBasket(Globals.userSession(), request.getParameter("basket-name")), list.size(), BasketActionType.REMOVEENTRY));
		}

		if (!Globals.session().contains(returnedObject))
			returnedObject = new ExperimentalProperty();

		return new WebModel(returnedObject).getModelAndView();

	}

	public boolean isBatchAction(String action)
	{
		return "addbasket".equals(action) 
				|| "removebasket".equals(action) 
				|| "addselect".equals(action) 
				|| "removeselect".equals(action)
				|| "publish".equals(action)
				|| ACTION_PUBLISH_FREE.equals(action);
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		return new WebModel().setList(Globals.getTaginationFilters(null)).getModelAndView();
	}

	private static Logger logger = LogManager.getLogger(EPBrowserController.class);

}

class Mutable<T>
{
	public T value;
}
