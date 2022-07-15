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

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.Hibernate;
import org.hibernate.HibernateException;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.business.WebFilters;
import qspr.entities.Article;
import qspr.entities.Basket;
import qspr.entities.BatchEditCompressedEP;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.MoleculeName;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.PropertyValue;
import qspr.entities.Unit;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.AccessChecker;
import qspr.util.CriteriaWrapper;
import cern.colt.Timer;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.mailer.Mailer;

@SuppressWarnings("unchecked")
@Controller
public class BatchEditController extends EPBrowserController 
{
	private static transient final Logger logger = LogManager.getLogger(BatchEditController.class);

	private static final int restartAfterNoRecords = 50;
	private static int VERBOSE = 0;

	public BatchEditController()
	{
		sessionRequired = true;
	}

	public ModelAndView show (HttpServletRequest req, HttpServletResponse res) throws Exception
	{
		return new WebModel().setTemplate("batchedit-browser").getModelAndView();	
	}

	/**
	 * get the filtered records and compress the (always)
	 * store the compression for the case when nothing was edited
	 */
	@Override
	public ModelAndView list(HttpServletRequest request, HttpServletResponse res) throws Exception
	{
		// get the selected records (these with a checked check box)
		Basket selectionBasket = Globals.selectionBasket(assertParam("trash"));
		if (selectionBasket.getRowsSize() == 0) // including "-1"
			throw new UserFriendlyException("You did not select any records. Please, go back and select the records you want to modify");
		//		
		// get the other filters
		WebFilters filters = formFilters(request);
		if (filters.has("basket-select") && ! filters.get("basket-select").equals(selectionBasket.id.toString()) && ! filters.get("basket-select").equals("-1"))
		{
			filters.addFilterOverride("basket-select", filters.get("basket-select") + "," + selectionBasket.id, null);
			//throw new UserFriendlyException("Batch editing while selecting another basket is not yet supported, coming soon");
		}

		// add selected records to the filters
		filters.unset("basket-select");
		filters.addFilter("basket-select", selectionBasket.id.toString(), null);

		CriteriaWrapper c = createMainCriteria(filters, Globals.RIGHTS_FREELY_AVAILABLE, Globals.session().createCriteria(ExperimentalProperty.class), new Mutable<Boolean>());
		c.criteria.setProjection(Projections.id());

		List<Long> filteredIDs = c.criteria.list();
		// List<Long> filteredIDs = Globals.selectionBasket(false).getIds();

		Set<Long> filteredSelectionIdSet = new HashSet<Long>();
		filteredSelectionIdSet.addAll(filteredIDs);
		filteredSelectionIdSet.remove(-1L); // why to remove -1

		Globals.setSessionAttribute(SessionVariable.FILTERED_SELECTION_IDSET, filteredSelectionIdSet);

		Map<String, BatchEditCompressedEP> beCompressedEPs = epCompression();
		Globals.setSessionAttribute(SessionVariable.COMPRESSED_EP_MAP, beCompressedEPs);

		WebList newWebList = new WebList().loadFromList(new ArrayList<BatchEditCompressedEP>(beCompressedEPs.values()));


		BrowserModel browserModel = new BrowserModel().setFilters(filters);
		return browserModel.setObject(newWebList).getModelAndView();

	}

	public ModelAndView edit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		String key = request.getParameter("key");
		List<BatchEditCompressedEP> compressedEPs = new ArrayList<BatchEditCompressedEP>(); // list of compressed eps sent for editing, should contain only one compressed ep
		Map<String, BatchEditCompressedEP> compressed = (Map<String, BatchEditCompressedEP>) Globals.getSessionAttribute(SessionVariable.COMPRESSED_EP_MAP);
		compressedEPs.add(compressed.get(key));
		WebList newWebList = new WebList().loadFromList(compressedEPs);
		return new WebModel(newWebList).setTemplate("batchedit-edit").getModelAndView();
	}

	public static List<Long> getBatchEditIds(Long propertyId, Long articleId, Long unitId) throws HibernateException, NumberFormatException
	{
		Criteria batchCriteria = Globals.session().createCriteria(ExperimentalProperty.class);
		Set<Long> selectionList = new HashSet<Long>();
		selectionList.addAll(Globals.selectionBasket(false).getIds());//Globals.userSession().selectionList;
		selectionList.remove(-1L);
		//		Set<Long> selectionList = Globals.userSession().selectionList;
		//		selectionList.remove(-1L);
		if (selectionList == null || selectionList.size() == 0)
			batchCriteria.add(Restrictions.sqlRestriction("1 = 0"));
		else
			batchCriteria.add(Restrictions.in("id", selectionList));

		batchCriteria.add(Restrictions.eq("property.id", propertyId))
		.add(Restrictions.eq("article.id", articleId));
		if(unitId != null)
			batchCriteria.add(Restrictions.eq("unit.id", unitId));
		else
			batchCriteria.add(Restrictions.isNull("unit.id"));

		batchCriteria.setProjection(Projections.distinct(Projections.id()));
		List<Long> basketEn =  batchCriteria.list();
		return basketEn;
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		String action = request.getParameter("action");
		WebList list = new WebList();
		list.list = new ArrayList<Object>();

		// if nothing is selected skip the action and return
		if ( ! assertParam("apply_prop") && 
				! assertParam("apply_page") &&
				! assertParam("apply_line") &&
				! assertParam("apply_table") &&
				! assertParam("apply_evi") &&
				! assertParam("apply_cond") &&
				! assertParam("apply_article") &&
				! assertParam("apply_public")
				)
		{
			return new WebModel(list).getModelAndView();
		}

		if (action.equals("edit"))
		{
			// get a list of experimental property IDs. These are the records to be updated
			Long propertyId = getLongParam("property");
			Long articleID = getLongParam("article");
			Long unitId = getLongParam("unit") != null ? getLongParam("unit") : 9L;

			Property property = Property.getById(propertyId);

			if(property.introducer == null)
				throw new UserFriendlyException("The property \"" + property.getName() + 
						"\" does not have yet introducer and description. Please, edit and save it before you can edit the records.");

			// epTemplate is the template ep after editing
			// modify the stored compressed ep and use it as template
			String usedKey = propertyId + "_" + articleID + "_" + unitId;
			//BatchEditCompressedEP becEP = getBatchEditCompressedEP(request, usedKey);

			BatchEditCompressedEP becEP = ((Map<String, BatchEditCompressedEP>) Globals.getSessionAttribute(SessionVariable.COMPRESSED_EP_MAP)).get(usedKey);
			List<Long> expPropIdList = new ArrayList<Long>(becEP.compresseEpIDs);

			ExperimentalProperty epTemplate = becEP.ep;
			updateTemplateEP(becEP, request); 				// modify compressed ep and use it as template
			Globals.session().evict(epTemplate);
			List<ExperimentalProperty> errors = new ArrayList<ExperimentalProperty>();
			List<Long> error_ep_ids = new ArrayList<Long>();
			int errorcount = 0;
			int successful = 0;

			Timer globalTimer = new Timer(); globalTimer.start();
			Timer localTimer = new Timer();
			while (expPropIdList.size() > 0)
			{
				localTimer.start();

				// get a sub list of IDs, after which the session will be restarted
				int highIndex = Math.min(restartAfterNoRecords, expPropIdList.size());
				List<Long> slice = expPropIdList.subList(0, highIndex);
				List<ExperimentalProperty> expPropList = 
						Globals.session().createCriteria(ExperimentalProperty.class)
						.add(Restrictions. in("id", slice))
						.list();
				slice.clear(); //Delete ids from original list

				// update the records
				for (ExperimentalProperty experimentalProperty : expPropList) 
				{
					AccessChecker.requestModificationPermission(experimentalProperty);
					if (assertParam("apply_public") && "true".equals(request.getParameter("apply_public"))){
						if(!property.isPublished()){

							Mailer.notifyAdmins("request to publish data for " + property.getName(), 
									" User: " + Globals.userSession().user + " for record: " + experimentalProperty.id);

							throw new UserFriendlyException("The property \"" + property.getName() + 
									"\" is hidden or not approved. Please, make it public and request administrator to approve it before you can publish the records.");
						}
						experimentalProperty.rights = Globals.RIGHTS_FREELY_AVAILABLE;
					}

					//set the edited property
					if (assertParam("apply_prop") && "true".equals(request.getParameter("apply_prop"))) 
					{
						if (epTemplate.property != null)
							experimentalProperty.property = epTemplate.property;
						if (epTemplate.unit != null && ! epTemplate.unit.multi)
							experimentalProperty.unit = epTemplate.unit;
						if (epTemplate.option != null && ! epTemplate.option.multi) {
							experimentalProperty.option = epTemplate.option;
							//switched of for the moment! experimentalProperty.value = epTemplate.value;
						}

						// predicate, secondValue
					}

					// set article
					if (assertParam("apply_article") && "true".equals(request.getParameter("apply_article"))) 
					{
						if (epTemplate.article != null)
							experimentalProperty.article = epTemplate.article;
					}

					// set article line, page, table
					if (assertParam("apply_line") && "true".equals(request.getParameter("apply_line"))) 
					{
						if (epTemplate.artLineNum != null) // delete this if, if empty field should remove all line numbers
							experimentalProperty.artLineNum = epTemplate.artLineNum;
					}
					if (assertParam("apply_page") && "true".equals(request.getParameter("apply_page"))) 
					{
						if (epTemplate.artPageNum != null) // delete this if, if empty field should remove all page numbers
							experimentalProperty.artPageNum = epTemplate.artPageNum;
					}
					if (assertParam("apply_table") && "true".equals(request.getParameter("apply_table"))) 
					{
						if (epTemplate.artTableNum != null) // delete this if, if empty field should remove all table numbers
							experimentalProperty.artTableNum = epTemplate.artTableNum;
					}

					// set status and evidence
					if (assertParam("apply_evi") && "true".equals(request.getParameter("apply_evi"))) 
					{
						if (null != becEP.ep_evidence )
						{	
							//experimentalProperty.ep_status = epTemplate.ep_status;

							switch (becEP.ep_evidence)
							{
							case -1: // nothing specified -> evidence should stay // do not change the original
							{
								break;
							}
							case 0: // no evidence specified
							{
								experimentalProperty.connectedProperty = null; 
								experimentalProperty.ep_status = null;
								break;
							}
							case 1: // measured in this article
							{
								experimentalProperty.connectedProperty = experimentalProperty; 
								experimentalProperty.ep_status = null;
								break;
							}
							case 2: // measured in this article (to be verified)
							{
								experimentalProperty.connectedProperty = experimentalProperty; 
								experimentalProperty.ep_status = ExperimentalProperty.STATUS_TOVERIFY;
								break;
							}
							case 3: // error in this record
							{
								experimentalProperty.ep_status = ExperimentalProperty.STATUS_ERROR;
								break;
							}
							case 4: // invalid record
							{
								experimentalProperty.ep_status = ExperimentalProperty.STATUS_INVALID;
								break;
							}
							}
						}
					}

					// conditions
					if (assertParam("apply_cond") && "true".equals(request.getParameter("apply_cond"))) 
					{
						if (epTemplate.conditions != null) {
							if(epTemplate.conditions.values.size() >0)
								experimentalProperty.conditions =
								updateConditionSet(epTemplate.conditions, experimentalProperty.conditions, request);
							else
								experimentalProperty.conditions = null;
						}
					}

					// Ownership and rights		
					if (Globals.userSession().user != null)
					{
						experimentalProperty.owner = Globals.userSession().user;

						if (assertParam("rights"))
						{
							experimentalProperty.rights = Integer.valueOf(request.getParameter("rights"));
							Globals.userSession().user.defaultRights = experimentalProperty.rights;
						}
					}

					// check for correctness, especially obligatory conditions

					experimentalProperty.firstEntry = null;
					try
					{
						experimentalProperty.updateHash();
						experimentalProperty.resolveConflictsAndSave();
						successful++;

					} catch (UserFriendlyException e) {
						throw e;
					} catch (Exception e)
					{
						experimentalProperty.ep_status = ExperimentalProperty.STATUS_ERROR;
						experimentalProperty.errorComment = e.getMessage();

						errorcount++;
						errors.add(experimentalProperty);
						error_ep_ids.add(experimentalProperty.id);
					}

					Globals.session().flush();

				}

				if (VERBOSE > 0) logger.info("Restarting transaction...");
				Globals.restartAllTransactions(true);

				if (VERBOSE > 0) logger.info("Per transaction: "+localTimer.seconds()+" sec");
				localTimer.reset();
				if (VERBOSE > 0) logger.info("Total: "+globalTimer.minutes()+" min");

			}
			logger.info(successful + " were saved susuccessfully");
			logger.info(errorcount + " could not be updated");

		}

		//TODO re-compress here and show the updated list

		return new WebModel().getModelAndView();
	}

	/**
	 * this template is build to apply the changes to all records, so it's formed after clicking the save button
	 * @param becEP 
	 * 
	 * @param request
	 * @param expPropIdList
	 * @return
	 * @throws Exception 
	 */
	private void updateTemplateEP(BatchEditCompressedEP becEP, HttpServletRequest request) throws Exception
	{
		String[] newkey = becEP.key.split("_"); 			 // get the current cache key and update it if necessary
		//get the new property and unit
		if (assertParam("apply_prop") && request.getParameter("apply_prop").equals("true"))
		{
			if (assertParam("newproperty")) 
			{
				Long propertyID = Long.valueOf(request.getParameter("newproperty"));
				Property property = (Property) Globals.session().get(Property.class, propertyID);
				Hibernate.initialize(property.options);
				Hibernate.initialize(property.obligatoryConditions);
				Globals.session().evict(property);
				becEP.ep.property = property;
				newkey[0] = "" + propertyID;
			}
			if (assertParam("newunit"))
			{
				Long unitID = Long.valueOf(request.getParameter("newunit"));
				if (unitID.longValue() != 0) { // unitID == 0 --> do not change
					if (unitID.longValue() < 0) // unit category has changed --> select default unit	
						unitID = becEP.ep.property.defaultUnit.id;
					becEP.ep.unit = (Unit) Globals.session().get(Unit.class, unitID);
					newkey[2] = "" + unitID;
				}
				// else do not change event

			}
			if (assertParam("newoption"))
			{
				Long optionID = Long.valueOf(request.getParameter("newoption"));
				if (optionID > 0) // -1 --> do not change
				{
					becEP.ep.option = (PropertyOption) Globals.session().get(PropertyOption.class, optionID);
					//switched off for now! becEP.ep.value = 0d; // reset the value to zero when the type of property changes
				}
			}
		} 

		// get the new article
		if (assertParam("apply_article") && request.getParameter("apply_article").equals("true"))
			if (assertParam("newarticle"))
			{
				Long articleID = Long.valueOf(request.getParameter("newarticle"));
				Article article = (Article) Globals.session().get(Article.class, articleID);
				Hibernate.initialize(article.pdfs);
				Hibernate.initialize(article.authors);
				Globals.session().evict(article);
				becEP.ep.article = article;
				newkey[1] = "" + articleID;
			}

		//get the new line, table and page number
		if (assertParam("apply_page") && request.getParameter("apply_page").equals("true") && assertParam("newpage"))
			becEP.ep.artPageNum = getIntParam("newpage");
		if (assertParam("apply_line") && request.getParameter("apply_line").equals("true") && assertParam("newline"))
			becEP.ep.artLineNum = Integer.valueOf(request.getParameter("newline"));
		if (assertParam("apply_table") && request.getParameter("apply_table").equals("true") && assertParam("newtable"))
			becEP.ep.artTableNum = request.getParameter("newtable");

		//get the edited evidence
		if (request.getParameter("apply_evi") != null && request.getParameter("apply_evi").equals("true"))
			if (assertParam("newevidence"))
				becEP.ep_evidence = getIntParam("newevidence"); 

		// update the conditionSet of the compressed ep to use it as template
		if (assertParam("apply_cond") && request.getParameter("apply_cond").equals("true"))
		{
			ConditionSet templateCSfromEditor = createTemplateConditionSetFromEditor(request);
			becEP.ep.conditions = updateTemplateConditionSet(templateCSfromEditor, becEP.ep.conditions, request);
			Collections.sort(becEP.ep.conditions.values, PropertyValue.condNameComp);
			//Globals.session().evict(becEP.ep.conditions);
		}

		String newKey = newkey[0] + "_" + newkey[1] + "_" + newkey[2];
		becEP.key = newKey;
		Globals.session().evict(becEP.ep); 
	}

	private ConditionSet createTemplateConditionSetFromEditor(HttpServletRequest request)
	{
		// create condition set
		String[] condition_ids = request.getParameterValues("cond-id");
		String[] old_condition_ids = request.getParameterValues("old-cond-id");
		String[] cond_values = request.getParameterValues("cond-value");
		String[] cond_unit_ids = request.getParameterValues("cond-unit");
		String[] cond_option_ids = request.getParameterValues("cond-option");
		String[] conditionTextValues = request.getParameterValues("cond-text-value");

		ConditionSet newConditionSet = new ConditionSet(); 		// new template ConditionSet

		if (condition_ids != null)
		{
			/////////////////////////////
			// Create new condition set
			for (int i = 0; i < condition_ids.length; i++) 
			{
				PropertyValue newConditionValue = new PropertyValue();

				newConditionValue.property = (Property) Globals.session().get(Property.class, Long.valueOf(condition_ids[i]));
				if (newConditionValue.property != null)
					Hibernate.initialize(newConditionValue.property.options);
				else
					throw new RuntimeException(condition_ids[i] + " is not vaid for condition_ids[i]!\nwhere does the -1 come from.");

				if (newConditionValue.property.isNumeric()) // quantitative condition
				{
					// value null means later: do not update
					// might be null, is used for discrimination later
					newConditionValue.value = ! cond_values[i].equals("") ? Double.valueOf(cond_values[i]) : null;
					newConditionValue.unit = (Unit) Globals.session().get(Unit.class, Long.valueOf(cond_unit_ids[i]));
				}
				else if (newConditionValue.property.isQualitative()) // qualitative Condition
				{
					// id -1 means later: do not update
					PropertyOption option = (PropertyOption) Globals.session().get(PropertyOption.class, Long.valueOf(cond_option_ids[i]));
					if (null != option)
					{
						newConditionValue.property = option.property;
						newConditionValue.option = option;
					} // otherwise property and option stays null 

					newConditionValue.unit = newConditionValue.property.defaultUnit;
				}
				else if (newConditionValue.property.isTextual())
				{
					// Textual condition
					if ( ! conditionTextValues[i].equals(""))
						newConditionValue.textualValue = conditionTextValues[i];
					// else its null
				}
				else
					throw new RuntimeException("Not implemeted yet");

				if (newConditionValue != null)
				{
					Hibernate.initialize(newConditionValue.property.unitCategory.units); // maybe eviction is necessary as well, so far not. ask Robert

					// no change of property if newConditionValue.property.id == newConditionValue.old_id 
					// else property has changed
					newConditionValue.old_id = Long.parseLong(old_condition_ids[i]);

					newConditionSet.values.add(newConditionValue);
				}
			}
		}
		return newConditionSet;
	}

	private PropertyValue getPropertyValue(Long prop_id, ConditionSet cs)
	{
		for (PropertyValue pv : cs.values)
			if (pv.property.id.equals(prop_id))
				return pv;

		return null;
	}

	private ConditionSet updateTemplateConditionSet(ConditionSet templateFromEditor, ConditionSet toUpdateConditionSet, HttpServletRequest request)
	{
		for (PropertyValue fromEditorCV : templateFromEditor.values)
		{
			PropertyValue toUpdateConditionPropertyValue = getPropertyValue(fromEditorCV.old_id, toUpdateConditionSet);

			if (toUpdateConditionPropertyValue != null) 				// update existing conditions
			{
				fromEditorCV.multi = toUpdateConditionPropertyValue.multi;
				toUpdateConditionSet.values.remove(toUpdateConditionPropertyValue);
				toUpdateConditionSet.values.add(fromEditorCV);
			}
			else														// set it to the result from editor (create non-existing conditions)
			{
				toUpdateConditionPropertyValue = fromEditorCV;
				toUpdateConditionSet.values.add(toUpdateConditionPropertyValue);  
			}
		}

		// delete conditions that have been removed
		List<PropertyValue> toRemove = new ArrayList<PropertyValue>();
		for (PropertyValue toUpdateCondition : toUpdateConditionSet.values)
		{
			// delete removed conditions
			if ( templateFromEditor.getPropertyValue(toUpdateCondition.property.id) == null )
			{
				toRemove.add(toUpdateCondition);
			}
		}
		toUpdateConditionSet.values.removeAll(toRemove);

		return toUpdateConditionSet;
	}

	/**
	 * This is the update for the real records. works slightly different than the update of the template.
	 * 
	 * @param templateSet: the template condition set, contains the changes to be applied
	 * @param expConditions: condition set to be updated
	 * @param request: 
	 * @throws Exception
	 */
	private ConditionSet updateConditionSet(ConditionSet templateSet, ConditionSet existingConditionSet, HttpServletRequest request) throws Exception
	{
		// new set, will replace expConditions set in the end
		ConditionSet updateSet = new ConditionSet();				
		//		if (existingConditionSet == null)
		//			existingConditionSet = new ConditionSet(); 				// to be able to loop over it

		// loop over template conditions		
		for (PropertyValue templateCV : templateSet.values) 		// template conditions
		{	
			String apply_id = "apply_cond_one_" + templateCV.property.id;
			boolean applyChange = assertParam(apply_id) && request.getParameter(apply_id).equals("on");

			// get the correct condition that was modified
			PropertyValue existingCV;

			if (existingConditionSet == null)
				existingCV = null;
			else
				existingCV = existingConditionSet.getPropertyValue(templateCV.old_id);

			if ( ! applyChange)
			{
				// add the existing condition untouched to the new set
				if (existingCV != null)
					updateSet.values.add(PropertyValue.clone(existingCV));
				continue;
			}

			// the template property is set in any case (same or new property)
			PropertyValue newCV = new PropertyValue();
			newCV.property = (Property) Globals.session().get(Property.class, Long.valueOf(templateCV.property.id));

			// modify existing condition or create it else
			if (existingCV != null)
			{
				// check type of existing condition  
				if (existingCV.property.isNumeric()) 			
				{	
					// check type of condition from the editor
					if (templateCV.property.isNumeric()) /// new CV or old CV is also textual, check if text should be updated
					{
						// change value or keep existing one
						if (null != templateCV.value)				 	
							newCV.value = templateCV.value;			// value field is not empty -> change			
						else 										
							newCV.value = existingCV.value;			// value field is empty (multi or by hand) -> no change

						// change unit or use existing one
						if (null != templateCV.unit) // change to new unit
							newCV.unit = (Unit) Globals.session().get(Unit.class, Long.valueOf(templateCV.unit.id));
						else
							newCV.unit = existingCV.unit; 			// use existing unit

						assert newCV.unit != null;
					}
					else if (templateCV.property.isQualitative())	// condition type has changed, apply defaults where possible
					{
						if (null != templateCV.option)
						{	
							// a new qualitative condition was selected, use the selected option
							if (templateCV.option.id != null  &&  templateCV.option.id != -1)
								newCV.option = (PropertyOption) Globals.session().get(PropertyOption.class, Long.valueOf(templateCV.option.id));
							assert newCV.option != null;
						} 
					}
					else if (templateCV.property.isTextual())
					{
						// a new textual condition was selected
						if (null != templateCV.textualValue)			// text field is not empty -> change
							newCV.textualValue = templateCV.textualValue;
						else // check in editor if filled
							newCV.textualValue = "no text specified";
					}
				}
				else if (existingCV.property.isQualitative())
				{
					if (templateCV.property.isQualitative()) /// new CV or old CV is also textual, check if text should be updated
					{
						if (null != templateCV.option) 	// change to new option
						{
							// a new qualitative condition was selected, use the selected option
							if (templateCV.option.id != null  &&  templateCV.option.id != -1)
								newCV.option = (PropertyOption) Globals.session().get(PropertyOption.class, Long.valueOf(templateCV.option.id));
							assert newCV.option != null;
						}
						else 										// keep the old option
						{
							//newCV = new PropertyValue(existingCV.option);
							newCV.option = (PropertyOption) Globals.session().get(PropertyOption.class, Long.valueOf(existingCV.option.id));
							assert newCV.option != null;
						}
					}
					else if (templateCV.property.isNumeric())
					{
						if (null != templateCV.value) 					// create condition with given value
							newCV.value = templateCV.value;
						// else value is null, this conditionValue will be filtered out in cv.isValid()
						// and therefore not created

						if (null != templateCV.unit) 					// put given unit
							newCV.unit = (Unit) Globals.session().get(Unit.class, Long.valueOf(templateCV.unit.id));
						else 										// put default unit
							newCV.unit = newCV.property.defaultUnit;

						assert newCV.unit != null;
					}
					else if (templateCV.property.isTextual())
					{
						// a new textual condition was selected
						if (null != templateCV.textualValue)			// text field is not empty -> change
							newCV.textualValue = templateCV.textualValue;
						else // check in editor if filled
							newCV.textualValue = "no text specified";
					}
				}
				else if (existingCV.property.isTextual())			// change to new textual condition
				{
					if (templateCV.property.isTextual()) /// new CV or old CV is also textual, check if text should be updated
					{
						if (null != templateCV.textualValue)			// text field is not empty -> change
							newCV.textualValue = templateCV.textualValue;
						else ///	soll zu werden new PropertyValue(tempPV.property, existing_cv.textualValue);									// text field is empty (multi or by hand) -> no change
							newCV.textualValue = existingCV.textualValue;	
					}
					else if (templateCV.property.isNumeric())
					{
						if (null != templateCV.value) 					// create condition with given value
							newCV.value = templateCV.value;
						// it should be checked by UI that the field is filled
						// else value is null, this conditionValue will be filtered out in cv.isValid()
						// and therefore not created

						if (null != templateCV.unit) 					// put given unit
							newCV.unit = (Unit) Globals.session().get(Unit.class, Long.valueOf(templateCV.unit.id));
						else 										// put default unit
							newCV.unit = newCV.property.defaultUnit;

						assert newCV.unit != null;
					}
					else if (templateCV.property.isQualitative())
					{
						// a new qualitative condition was selected, use the selected option
						if (null != templateCV.option)
						{
							if (templateCV.option.id != null  &&  templateCV.option.id != -1)
								newCV.option = (PropertyOption) Globals.session().get(PropertyOption.class, Long.valueOf(templateCV.option.id));
							assert newCV.option != null;
						}
					}
				}
				else {
					throw new RuntimeException("Not implemeted yet");
				}
			}
			else // create this condition
			{
				if (templateCV.property.isNumeric())
				{
					if (null != templateCV.value) 					// create condition with given value
						newCV.value = templateCV.value;
					// else value is null, this conditionValue will be filtered out in cv.isValid()
					// and therefore not created

					if (null != templateCV.unit) 					// put given unit
						newCV.unit = (Unit) Globals.session().get(Unit.class, Long.valueOf(templateCV.unit.id));
					else 										// put default unit
						newCV.unit = newCV.property.defaultUnit;

					assert newCV.unit != null;
				}
				else if (templateCV.property.isQualitative())
				{
					if (null != templateCV.option)
					{
						if (templateCV.option.id != null  &&  templateCV.option.id != -1)
							newCV.option = (PropertyOption) Globals.session().get(PropertyOption.class, Long.valueOf(templateCV.option.id));
						assert newCV.option != null;
					} 
				}
				else if (templateCV.property.isTextual())			// change to new textual condition
				{
					if (null != templateCV.textualValue)			// text field is not empty -> change
						newCV.textualValue = templateCV.textualValue;
					else // check in editor if filled
						newCV.textualValue = "no text specified";
				}
				else 
					throw new RuntimeException("Not implemeted yet");
			}


			if (newCV.isValid()) 								// check for complete condition values
				updateSet.values.add(newCV);
		}

		return updateSet.get();

	}

	private Map<String, BatchEditCompressedEP> epCompression() throws Exception
	{
		Set<Long> selectionSet = (Set<Long>) Globals.getSessionAttribute(SessionVariable.FILTERED_SELECTION_IDSET);
		Map<String, List<Long>> groups = groupBy_PropId_ArtId_UnitId(selectionSet);

		List<String> orderedKeys = new ArrayList<String>(groups.keySet());

		if (VERBOSE > 0) logger.info("****");
		if (VERBOSE > 0) logger.info(orderedKeys.toString());
		if (VERBOSE > 0) logger.info("****");

		Collections.sort(orderedKeys);

		Map<String, BatchEditCompressedEP> beCompressedEPs = new HashMap<String, BatchEditCompressedEP>();
		for (String key : orderedKeys) 
		{
			List<Long> props = groups.get(key);
			beCompressedEPs.put(key, compressGroupToEP(key, props));
		}
		return beCompressedEPs;
	}

	private Map<String, List<Long>> groupBy_PropId_ArtId_UnitId(Set<Long> expPropIds) throws Exception
	{
		Map<String, List<Long>> newMap = new HashMap<String, List<Long>>();
		String key = null;

		Timer t = new Timer(); t.reset(); t.start();

		List<ExperimentalProperty> ep_block;

		if (expPropIds.size() > 0)
			ep_block = Globals.session().createCriteria(ExperimentalProperty.class)
			.add(Restrictions.in("id", expPropIds))
			.list();
		else
			ep_block = new ArrayList<ExperimentalProperty>();

		//group selected set according to similar property, article and unit.
		for (ExperimentalProperty ep : ep_block) 
		{
			key = ep.property.id + "_" + ep.article.id + "_" + ((ep.unit != null) ? ep.unit.id : "9");
			if (VERBOSE > 0) logger.info(key);
			if (newMap.containsKey(key)) {
				newMap.get(key).add(ep.id);
			} else {
				List<Long> idList = new ArrayList<Long>();
				idList.add(ep.id);
				newMap.put(key, idList);
			}
		}
		t.stop(); if (VERBOSE > 0) logger.info("grouping set " + t.seconds()); 
		return newMap;
	}

	private int idCounter = 0;

	/**
	 * The compressed group is the compression of the records before they are edited in the editor
	 * This has to be done only once and then stored in the session until the selection of records changes
	 * 
	 * @param compressEpGroups
	 * @param key
	 * @param expPropIds
	 * @return BatchEditCompressedEP
	 */
	private BatchEditCompressedEP compressGroupToEP(String key, List<Long> expPropIds) {

		BatchEditCompressedEP becEP = new BatchEditCompressedEP();
		becEP.id = Long.valueOf(idCounter++);
		becEP.key = key;
		becEP.imageId = new ArrayList<Long>();

		becEP.ep = new ExperimentalProperty();
		becEP.ep.conditions = new ConditionSet();
		becEP.ep.moleculenames = new ArrayList<MoleculeName>();

		Map<Long, PropertyValue> conditionCompressor = new HashMap<Long, PropertyValue>(); 
		// get data from first experimentalProperty in list
		if (expPropIds.size() > 0)
		{

			ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, expPropIds.remove(0));
			becEP.compresseEpIDs.add(ep.id);

			becEP.imageId.add(ep.molecule.id);

			becEP.ep.article = ep.article;
			becEP.ep.property = ep.property;

			Hibernate.initialize(becEP.ep.property.obligatoryConditions);
			Hibernate.initialize(becEP.ep.article.pdfs);
			Hibernate.initialize(becEP.ep.article.authors);

			becEP.ep.value = ep.value;
			becEP.ep.unit = ep.unit;
			becEP.ep.option = ep.option;

			becEP.ep.artLineNum = ep.artLineNum;
			becEP.ep.artPageNum = ep.artPageNum;
			becEP.ep.artParagraph = ep.artParagraph;
			becEP.ep.artTableNum = ep.artTableNum;

			becEP.ep.ep_status = ep.ep_status;
			becEP.ep_evidence = determineEvidence(ep);

			becEP.ep.rights = ep.rights;
			if(ep.conditions != null)
			{
				for (PropertyValue pv : ep.conditions.values)
				{
					conditionCompressor.put(pv.property.id, createNewProperty(pv));
				}	
			}
		}

		if (expPropIds.size() > 0)
		{
			// get data from further experimentalProperty in list
			List<ExperimentalProperty> ep_block = Globals.session().createCriteria(ExperimentalProperty.class).add(Restrictions.in("id", expPropIds)).list();
			for (ExperimentalProperty ep : ep_block) // maybe it is faster here 
			{
				becEP.compresseEpIDs.add(ep.id);
				becEP.imageId.add(ep.molecule.id);

				becEP.ep.value = (becEP.ep.value == ep.value) ? ep.value : -0.00d;

				if ( becEP.ep.option != null && ! becEP.ep.option.name.equals(ep.option.name)) 
					becEP.ep.option.multi = true;

				becEP.ep.artLineNum = (null != becEP.ep.artLineNum && null != ep.artLineNum 
						&& becEP.ep.artLineNum.compareTo(ep.artLineNum) == 0) 
						? ep.artLineNum : null;
				becEP.ep.artPageNum = (null != becEP.ep.artPageNum && null != ep.artPageNum 
						&& becEP.ep.artPageNum.compareTo(ep.artPageNum) == 0)
						? ep.artPageNum : null;
				becEP.ep.artParagraph = (null != becEP.ep.artParagraph && null != ep.artParagraph 
						&& becEP.ep.artParagraph.compareTo(ep.artParagraph) == 0) 
						? ep.artParagraph : null;
				becEP.ep.artTableNum = (null != becEP.ep.artTableNum && null != ep.artTableNum 
						&& becEP.ep.artTableNum.compareTo(ep.artTableNum) == 0)
						? ep.artTableNum : null;

				Integer ep_evidence = determineEvidence(ep);
				becEP.ep_evidence = (null != becEP.ep_evidence && null != ep_evidence && 
						becEP.ep_evidence.compareTo(ep_evidence) == 0)
						? ep_evidence : null;

				becEP.ep.rights = (becEP.ep.rights == ep.rights) ? ep.rights : Globals.RIGHTS_FREELY_AVAILABLE;

				if (ep.conditions != null)
				{
					for (PropertyValue expCv : ep.conditions.values) 
					{
						Long propertyKey = expCv.property.id;
						PropertyValue pv = conditionCompressor.get(propertyKey);
						if (null == pv) 			//new condition value
						{
							conditionCompressor.put(expCv.property.id, createNewProperty(expCv));
						}
						else 						// condition value has already been there
						{
							if (expCv.property.isQualitative())
							{
								if ( ! pv.option.id.equals(expCv.option.id))
									pv.multi = true;
							}
							else if (expCv.property.isNumeric())
							{
								if ( ! pv.value.equals(expCv.value))
									pv.multi = true;

								Unit cond_unit = (Unit) Globals.session().get(Unit.class, Long.valueOf(expCv.unit.id));
								if (cond_unit.id != pv.unit.id)
									pv.unit.multi = true;
							}
							else if (expCv.property.isTextual())
							{
								if ( ! expCv.textualValue.equals(pv.textualValue))
									pv.multi = true;
							}
							else
								throw new RuntimeException("Not implemeted yet");
						}
					}
				}
			}
		}

		// prepare a new condition set that contains all conditions that appear in any of the selected records
		// if they have different values, "---" is displayed
		ConditionSet cs = new ConditionSet();

		cs.values.addAll(conditionCompressor.values());
		Collections.sort(cs.values, PropertyValue.condNameComp);

		if (VERBOSE > 0)
			for (PropertyValue l : cs.values)
				logger.info(l.property.unitCategory.units.size() + " - " + l.toString());

		//request().getSession().setAttribute(COMPRESSED_CONDITION_SET, cs);		
		becEP.ep.conditions = cs;
		Globals.session().evict(becEP.ep);
		return becEP;
	}

	private PropertyValue createNewProperty(PropertyValue pv)
	{
		PropertyValue newPropertyValue = null;
		switch (pv.property.type)
		{
		case Property.TYPE_NUMERIC:
			Unit cond_unit = (Unit) Globals.session().get(Unit.class, Long.valueOf(pv.unit.id));
			newPropertyValue = new PropertyValue(pv.property, pv.value, cond_unit);
			break;
		case Property.TYPE_QUALITATIVE:
			PropertyOption option = (PropertyOption) Globals.session().get(PropertyOption.class, Long.valueOf(pv.option.id));
			newPropertyValue = new PropertyValue(option);
			break;
		case Property.TYPE_TEXTUAL:
			newPropertyValue = new PropertyValue(pv.property, pv.textualValue);
			break;
		}
		newPropertyValue.id = pv.property.id;
		Hibernate.initialize(newPropertyValue.property.unitCategory.units);
		return newPropertyValue;
	}


	private Integer determineEvidence(ExperimentalProperty ep) {
		if (null != ep.ep_status)
		{	
			if ( ep.ep_status == ExperimentalProperty.STATUS_ERROR )
				return 3;  // error in this record
			else if (ep.ep_status == ExperimentalProperty.STATUS_INVALID)
				return 4;  // invalid record
			else if (ep.ep_status == ExperimentalProperty.STATUS_TOVERIFY)
				return 2;  // measured in this article (to be verified)

		}
		else
		{
			if (null == ep.connectedProperty && ep.ep_status == null)
				return 0; // no evidence specified
			else if	(ep.id.longValue() == ep.connectedProperty.id.longValue() && ep.ep_status == null)
				return 1; // measured in this article
			//			else if (ep.id.longValue() == ep.connectedProperty.id.longValue() && ep.ep_status == ExperimentalProperty.STATUS_TOVERIFY)
			//				return 2; // measured in this article (to be verified)
		}
		return null;
	}


}

