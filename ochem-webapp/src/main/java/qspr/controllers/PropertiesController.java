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

import java.io.BufferedReader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.SessionVariable;
import qspr.ThreadScope;
import qspr.business.WebFilters;
import qspr.dao.PropertyDAO;
import qspr.dao.PropertyDAOImpl;
import qspr.dao.Repository;
import qspr.entities.Alert;
import qspr.entities.BatchEditCompressedEP;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.PropertyOptionsFilter;
import qspr.entities.PropertyValue;
import qspr.frontend.BrowserModel;
import qspr.frontend.LabeledValue;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;

import com.eadmet.business.PaginationFilter;
import com.eadmet.business.PropertiesAction;
import com.eadmet.business.PropertiesFilter;
import com.eadmet.business.PropertiesFilter.ApprovalStatus;
import com.eadmet.business.PropertiesService;
import com.eadmet.exceptions.UserFriendlyException;

@Controller
public class PropertiesController extends BrowserWrapper
{

	private static Logger logger = LogManager.getLogger(PropertiesController.class);

	// Autowire?
	PropertiesService service = new PropertiesService();
	PropertyDAO dao = new PropertyDAOImpl();

	public PropertiesController()
	{
		sessionRequired = true;
	}

	public ModelAndView autocomplete(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		PropertiesFilter filter = getPropertiesFilter(request);
		PaginationFilter pager = getPaginationFilter(request);
		List<LabeledValue> propertyLabels = service.getLabels(filter, pager);
		return new WebModel().addObjects(propertyLabels).getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		String action = request.getParameter("action");

		if (action == null)
			return new WebModel().getModelAndView();

		Property p = null;

		if (request.getParameter("id") != null)
			p = dao.getPropertyById(Long.valueOf(request.getParameter("id")));

		if ("delete".equals(action))
			service.delete(p);
		else if ("approve".equals(action))
			service.approve(p);
		else if ("unapprove".equals(action))
			service.unapprove(p);
		else if ("publish".equals(action))
			service.publish(p);
		else if ("removechild".equals(action))
			service.removechild(p);
		else if ("addchild".equals(action))
		{
			Property parent = dao.getPropertyById(Long.valueOf(request.getParameter("parent")));
			p = dao.getPropertyById(Long.valueOf(request.getParameter("child")));
			service.addchild(parent, p);
		} 		
		else if ("edit".equals(action))
		{
			PropertiesAction propertiesAction = getPropertiesAction(request);
			p = service.edit(propertiesAction);
		} else
			throw new UserFriendlyException("Unknown action "+action+" in properties controller");

		return new WebModel(p).getModelAndView();
	}

	public ModelAndView addmany(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		if (Globals.userSession().user == null)
			return redirect("user/pleaseregister.do");

		if (request.getParameter("ids") != null) {
			String s = request.getParameter("ids");

			BufferedReader br = new BufferedReader(new StringReader(s));
			List<String> names = new ArrayList<String>();
			String line;
			while ((line = br.readLine()) != null) {
				names.add(line.trim());
			}
			br.close();
			String options[] = new String[]{"high","low"};

			boolean newP = false;

			for(String name:names){
				Property prop = Repository.property.getProperty(name, false);
				if(prop == null) {
					prop = new Property();
					prop.setName(name);
					prop.type = Property.TYPE_QUALITATIVE;
					prop.defaultUnit = Repository.unit.getUnitById(289);
					prop.unitCategory = prop.defaultUnit.category;
					prop.isCondition = false;
					prop.owner = ThreadScope.get().userSession.user;
					newP = true;
				}
				for(String option:options)
					if(prop.getOptionByName(option) == null) {
						prop.options.add(Repository.option.getPropertyOptionByName(option, prop.id, true, false));
						newP = true;
					}

				if(newP) {
					System.out.println("add " + name);
					Globals.session().saveOrUpdate(prop);
				}
			}
		}

		return new WebModel().setRenderMode("popup")
				.setTemplate("property-multi")
				.getModelAndView();
	}

	public ModelAndView edit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		if (Globals.userSession().user == null)
			return redirect("user/pleaseregister.do");

		Property prop;
		if (request.getParameter("id") == null || request.getParameter("id").equals("-1"))
			prop = service.create(getPropertiesAction(request));
		else
			prop = dao.getPropertyById(Long.valueOf(request.getParameter("id")));

		//service.fillUsedConditions(prop); // We do it dynamically with AJAX for a better user experience

		List<Object> list = new ArrayList<Object>();
		list.addAll(service.getUnitCategories());
		list.addAll(prop.obligatoryConditions);

		Globals.setMarshallingOption(MarshallingOption.PROPERTY_RECORD_COUNT);
		Globals.setMarshallingOption(MarshallingOption.PROPERTY_OPTIONS);
		Globals.setMarshallingOption(MarshallingOption.PROPERTY_MODERATOR);	
		Globals.setMarshallingOption(MarshallingOption.PROPERTY_UNITS);

		return new WebModel(prop)
				.setList(list)
				.setRenderMode("popup")
				.setTemplate("property-edit")
				.getModelAndView();
	}

	public ModelAndView getUsedConditions(HttpServletRequest request, HttpServletResponse response) {
		Property prop = dao.getPropertyById(getLongParam("id"));
		service.fillUsedConditions(prop);
		return new WebModel(prop).getModelAndView();
	}

	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		PropertiesFilter filter = getPropertiesFilter(request);
		PaginationFilter pager = getPaginationFilter(request);

		WebFilters webFilters = formFilters(request); 

		List<Property> properties = service.get(filter, pager);

		WebList list = new WebList().loadFromPagedList(properties, pager.pageNum, pager.pageSize, pager.totalSize);

		if (assertParam("show-counts"))
			Globals.setMarshallingOption(MarshallingOption.PROPERTY_RECORD_COUNT);
		Globals.setMarshallingOption(MarshallingOption.PROPERTY_MODERATOR);

		return new BrowserModel().setFilters(webFilters).setObject(list).setTemplate("property-browser").getModelAndView();
	}

	public ModelAndView applytags(HttpServletRequest request, HttpServletResponse response)
	{
		// Undocumented feature, batched tags application by name mask
		// No user interface supplied
		// To speed it up...
		service.applyTags(request.getParameter("mask"), request.getParameter("tag"));
		return redirect("properties/show.do");
	}

	public ModelAndView listoptions(HttpServletRequest req, HttpServletResponse res)
	{
		PropertyOptionsFilter filter = getPropertyOptionsFilter(req);
		logger.info("Using filter: "+filter);
		List<PropertyOption> options = Property.getOptions(filter);
		for(PropertyOption o: options) // since we restart transaction just below...
			Globals.session().evict(o);
		Globals.restartAllTransactions(true);
		WebList list = new WebList().loadFromList(options);
		return new WebModel(list).getModelAndView();
	}

	public ModelAndView listpredicates(HttpServletRequest req, HttpServletResponse res)
	{
		return new WebModel(new WebList().loadFromList(service.getPredicates())).getModelAndView();
	}

	public ModelAndView listvalues(HttpServletRequest request, HttpServletResponse res)
	{
		Globals.setMarshallingOption(MarshallingOption.PROPERTY_PREDICATES);
		Globals.setMarshallingOption(MarshallingOption.PROPERTY_OPTIONS);
		Globals.setMarshallingOption(MarshallingOption.UNITCATEGORY_UNITS); // Load units of categories, they are required

		List<PropertyValue> values = new ArrayList<PropertyValue>();

		if (assertParam("id"))
		{
			ExperimentalProperty ep = (ExperimentalProperty)Globals.session().get(ExperimentalProperty.class, getLongParam("id"));
			if (ep != null)
				if (ep.conditions != null)
					values = ep.conditions.values;
		}

		if (assertParam("set-id"))
		{
			ConditionSet cset = (ConditionSet) Globals.session().get(ConditionSet.class, getLongParam("set-id"));
			if (cset != null)
				values = cset.values;
		}

		return new WebModel(new WebList().loadFromList(values)).getModelAndView();
		//		if(assertParam("proid"))
		//		{
		//			Long propertyId = Long.valueOf(request.getParameter("proid"));
		//			Long articleId = Long.valueOf(request.getParameter("artid"));
		//			Long unitId = Long.valueOf(request.getParameter("unitid"));
		//			List<Long> experimentalProperties_ids = BatchEditController.getBatchEditIds(propertyId, articleId, unitId);
		//			return new WebModel(new WebList().loadFromList(conditionsGroup(experimentalProperties_ids))).getModelAndView();
		//		}
		//		return null;
	}

	@SuppressWarnings("unchecked")
	public ModelAndView listvaluesbatch(HttpServletRequest request, HttpServletResponse res)
	{
		Globals.setMarshallingOption(MarshallingOption.UNITCATEGORY_UNITS);
		Globals.setMarshallingOption(MarshallingOption.PROPERTY_OPTIONS_FULL);

		String key = request.getParameter("key");
		Map<String, BatchEditCompressedEP> compressed = (Map<String, BatchEditCompressedEP>) Globals.getSessionAttribute(SessionVariable.COMPRESSED_EP_MAP);
		BatchEditCompressedEP beCompEP = compressed.get(key);

		ConditionSet cs = null;
		if (beCompEP != null)
			cs = beCompEP.ep.conditions;

		if (cs != null) {
			return new WebModel(new WebList().loadFromList(cs.values)).getModelAndView();
		} else
			return new WebModel().getModelAndView();
	}

	public ModelAndView saveoptions(HttpServletRequest req, HttpServletResponse res)
	{
		String[] ids = req.getParameterValues("co-id");
		String[] names = req.getParameterValues("co-name");
		Property property = dao.getPropertyById(getLongParam("id"));
		service.saveOptions(property, ids, names);
		return new WebModel(new Alert("ok")).getModelAndView();
	}

	//	private List<PropertyValue> conditionsGroup(List<Long> prop_ids) 
	//	{
	//		BatchEditCompressedEP md = new BatchEditCompressedEP();
	//		md.ep = new ExperimentalProperty();
	//		md.ep.conditions = new ConditionSet();
	//		int i = 0;
	//		
	//		
	//		final int step = 50;
	//		final int size = prop_ids.size();
	//		int start = 0;
	//		int end = (step < size) ? step : size;
	//		
	//		while (start < size)
	//		{
	//			List<Long> blockList = prop_ids.subList(start, end);
	//			start = end;
	//			end = (start + step < size) ? start + step : size;
	//			
	//			List<ExperimentalProperty> ep_block = Globals.session().createCriteria(ExperimentalProperty.class).add(Restrictions.in("id", blockList)).list();
	//			
	//			for (ExperimentalProperty ep : ep_block) {
	//				if(i==0)
	//				{	
	//					md.id = ep.id;
	//					if(ep.conditions != null)
	//					{
	//						md.ep.conditions.values.addAll(ep.conditions.values);
	//					}
	//					i++;
	//				}
	//				else
	//				{ 
	//					if(ep.conditions != null)
	//					{
	//						for (PropertyValue ecv : ep.conditions.values) {
	//							boolean isin = false;
	//							for (PropertyValue ucv : md.ep.conditions.values) {
	//								if (ecv.property.id == ucv.property.id){ // use comparator here
	//									isin = true;
	//									if (ucv.property.isQualitative())
	//									{
	//										if(!ucv.option.id.equals(ecv.option.id))
	//											ucv.multi = true;
	//									}
	//									else if (ucv.property.isNumeric())
	//									{
	//										if(!ucv.value.equals(ecv.value))
	//											ucv.multi = true;
	//									}
	//									else
	//										if (!ucv.textualValue.equals(ecv.textualValue))
	//											ucv.multi = true;
	//									break;
	//								}
	//							}
	//							if ( ! isin)
	//								md.ep.conditions.values.add(ecv);
	//						}
	//					}
	//				}
	//			}
	//			
	//		}
	//		return md.ep.conditions.values;
	//	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		WebModel wm =  new WebModel().setList(Globals.getTaginationFilters(Property.class)).setTemplate("property-browser");
		if (assertParam("condition"))
			wm.addParam("condition", "true");
		return wm.getModelAndView();
	}

	public static PropertiesFilter getPropertiesFilter(HttpServletRequest r) 
	{
		PropertiesFilter pf = new PropertiesFilter();

		if (r.getParameter("name") != null)
			pf.name = r.getParameter("name");

		if (r.getParameter("term") != null) // Autocomplete reverse compatibility
			pf.name = r.getParameter("term");

		if (r.getParameter("query") != null)
			pf.query = r.getParameter("query");

		pf.condition = (r.getParameter("condition") != null);

		if (r.getParameter("parent") != null)
			pf.parentId = Long.valueOf(r.getParameter("parent"));

		if (r.getParameter("directories") != null)
			if ("true".equals(r.getParameter("directories")))
				pf.directories = true;
			else if ("false".equals(r.getParameter("directories")))
				pf.directories = false;

		if (r.getParameter("basket") != null)
			pf.basketId = Long.valueOf(r.getParameter("basket"));

		if (r.getParameter("tag") != null)
			pf.tagId = Long.valueOf(r.getParameter("tag"));


		if (r.getParameter("approval-status") != null)
			if ("only-awaiting-approval".equals(r.getParameter("approval-status")))
				pf.approvalStatus = ApprovalStatus.AWAITING_ONLY;
			else if ("only-approved".equals(r.getParameter("approval-status")))
				pf.approvalStatus = ApprovalStatus.APPROVED_ONLY;

		if (r.getParameter("id") != null)
			pf.id = Long.valueOf(r.getParameter("id"));
		return pf;
	}

	public static PropertiesAction getPropertiesAction(HttpServletRequest r) 
	{
		PropertiesAction pa = new PropertiesAction();

		if (r.getParameter("id") != null)
			pa.id = Long.valueOf(r.getParameter("id"));

		if (r.getParameter("parent") != null)
			pa.parentId = Long.valueOf(r.getParameter("parent"));

		if (r.getParameter("name") != null)
			pa.name = r.getParameter("name");

		if (r.getParameter("is-directory") != null)
			pa.isDirectory = true;

		if (r.getParameter("type") != null)
			pa.propertyType = Integer.valueOf(r.getParameter("type"));	

		if (r.getParameter("condition") != null && r.getParameter("condition").equals("true"))
			pa.isCondition = true;

		if (r.getParameter("category") != null)
			pa.unitCategoryId = Long.valueOf(r.getParameter("category"));

		if (r.getParameter("confirmed") != null && r.getParameter("confirmed").equals("true"));
		pa.confirmed = true;

		pa.description = r.getParameter("description");
		pa.aliases = r.getParameter("aliases");

		if (r.getParameter("unit") != null)
			pa.defaultUnitId = Long.valueOf(r.getParameter("unit"));

		if (r.getParameter("public") != null)
			pa.isPublic = true;

		if (r.getParameter("bonusPointsWeight") != null)
			pa.bonusPointsWeight = Double.valueOf(r.getParameter("bonusPointsWeight"));

		if (r.getParameter("approve") != null)
			pa.isApproved = true;

		if (r.getParameter("condition-id") != null)
			for (String conditionId :  r.getParameterValues("condition-id"))
				pa.obligatoryConditionIds.add(Long.valueOf(conditionId));

		return pa;
	}

	public static PropertyOptionsFilter getPropertyOptionsFilter(HttpServletRequest r)
	{
		PropertyOptionsFilter pof = new PropertyOptionsFilter();

		if (r.getParameter("id") != null)
			pof.propertyId = Long.valueOf(r.getParameter("id"));

		if (r.getParameter("basket") != null)
			pof.basketId = Long.valueOf(r.getParameter("basket"));
		else
			if(Globals.getSessionAttribute(SessionVariable.BASKET_SELECT)!=null)
				pof.basketId = (Long)Globals.getSessionAttribute(SessionVariable.BASKET_SELECT);

		return pof;
	}
}


