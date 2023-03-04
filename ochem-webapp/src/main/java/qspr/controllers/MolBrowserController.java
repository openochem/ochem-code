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
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.business.WebFilters;
import qspr.dao.Repository;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Fragment;
import qspr.entities.Mapping1;
import qspr.entities.Molecule;
import qspr.entities.Property;
import qspr.entities.Tag;
import qspr.frontend.BrowserModel;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.CriteriaWrapper;

import com.eadmet.exceptions.UserFriendlyException;

@Controller
@SuppressWarnings("unchecked")
public class MolBrowserController extends BrowserWrapper 
{
	private static transient final Logger logger = LogManager.getLogger(MolBrowserController.class);

	public MolBrowserController()
	{
		sessionRequired = true;
	}

	protected CriteriaWrapper createMainCriteria(WebFilters filters, Criteria criteria)
	{
		CriteriaWrapper wrapper = new CriteriaWrapper(criteria);

		if (OCHEMConfiguration.hidePrivateMolecules)
			criteria.add(Restrictions.eq("visible", Boolean.TRUE));

		if (filters.has("id"))
		{
			// ID has been explicitly specified
			// No need to look for other filters

			String[] ids = ThreadScope.get().localRequest.getParameterValues("id");
			if (ids == null)
				ids = new String[]{filters.get("id")};

			List<Long> mp1List = new ArrayList<Long>();
			List<Integer> mp2List = new ArrayList<Integer>();
			if (ids.length == 1)
			{
				ids = ids[0].split(",");

				for (String id : ids)
				{
					id = id.trim().toUpperCase();
					if (id.startsWith("M"))
						mp2List.add(Integer.valueOf(id.substring(1)));
					else
						mp1List.add(Long.valueOf(id));
				}
			}
			else
			{
				for (String id : ids)
				{
					id = id.trim().toUpperCase();
					if (id.startsWith("M"))
						mp2List.add(Integer.valueOf(id.substring(1)));
					else
						mp1List.add(Long.valueOf(id));
				}
			}

			if (!mp1List.isEmpty())
				criteria.add(Restrictions.in("id", mp1List));
			if (!mp2List.isEmpty())
				criteria.createCriteria("mapping2").add(Restrictions.in("id", mp2List));
		}
		else
		{
			if (filters.has("with-records") || filters.has("name") || filters.has("type"))
			{
				wrapper.createAlias("molecules", "m");
			}
			if (filters.has("with-records"))
			{
				Criteria epCriteria = criteria.createCriteria("m.experimentalProperties");
				epCriteria.add(Restrictions.isNull("deleted"));
				ExperimentalProperty.addAccessRestrictions(epCriteria, Globals.RIGHTS_FREELY_AVAILABLE, Globals.userSession().user, false, true);
			}
			if(filters.has("name"))
			{
				String name = filters.get("name");
				if (name.matches("[A-Z]{14}"))
				{
					// InchiKey-1
					criteria.add(Restrictions.eq("inchi1", name));
				}
				else if (name.matches("[A-Z]{14}-[A-Z]{10}"))
				{
					// Inchikey-1+2
					String[] parts = name.split("-");
					criteria.createCriteria("mapping2")
					.add(Restrictions.eq("inchi1", parts[0]))
					.add(Restrictions.eq("inchi2", name));
				}
				else
				{
					criteria
					.createCriteria("m.experimentalProperties")
					.createCriteria("moleculenames")
					.add(Restrictions.like("name", "%"+name+"%"));
				}
			}
			if (filters.has("fragSearch"))
			{

				String fragSearch = filters.get("fragSearch");
				filters.unset("fragSearch");
				CriteriaWrapper cwMp1 = createMainCriteria(filters, Globals.session().createCriteria(Mapping1.class));

				//				if (!cwMp1.hasAlias("mp1"))
				//				{
				////					cwMp1.createAlias("molecule", "mol");
				//					cwMp1.createAlias("mapping1", "mp1");
				//				}
				cwMp1.criteria.setProjection(Projections.groupProperty("id"));
				List<Long> mp1ids = cwMp1.criteria.list();
				if (mp1ids.size() > 70000)
					throw new UserFriendlyException("<br/>Sorry, the fragment filter can't be used with more than 70000 molecules." +
							"<br/>(Currently, the fragment filter is applied to " + mp1ids.size() + " molecules)." +
							"<br/>Please narrow your search together with other filters."
							);

				logger.info("Frag search: sending " + mp1ids.size() + " molecules");
				filters.addFilter("fragSearch", fragSearch, null);

				// fragment search in molecule browser
				Molecule fragMol = Repository.molecule.getMolecule(Long.valueOf(filters.get("fragSearch")));
				// next line is important to fill temporary local database
				Fragment frag = Fragment.get(fragMol.mapping1.inchi1, mp1ids, filters.getHash());
				criteria.createCriteria("moleculeFragment")
				.add(Restrictions.eq("id", frag.id)); // fragment context
			}
			if(filters.has("frag_id"))
			{
				List<Long> converted_mp1_ids = new ArrayList<Long>();
				List<Integer> mapping1_ids = qspr.fragmententities.Mapping1.getMapping1IdListByFragmentID(Integer.valueOf(filters.get("frag_id")));
				for (Integer long1 : mapping1_ids)
				{
					converted_mp1_ids.add(Long.valueOf(long1.toString()));
				}
				criteria.add(Restrictions.in("id", converted_mp1_ids));

			}
			if(filters.has("type"))
			{
				if(filters.get("type").equals("property"))
				{
					List<Property> pls = Globals.session().createCriteria(Property.class)
							.createAlias("tags", "t")
							.add(Restrictions.eq("t.id", filters.getLong("tag"))).list();
					wrapper.createAlias("m.experimentalProperties", "ep")
					.add(Restrictions.in("ep.property", pls.toArray()));
				}else{
					criteria.createAlias("tags", "t")
					.add(Restrictions.eq("t.id", filters.getLong("tag")));
				}
			}
		}

		Globals.applyTaginationFilters(criteria, Mapping1.class, "mp1");

		criteria.addOrder(Order.desc("id"));
		return wrapper;

	}

	private void executeAction(Mapping1 molecule, HttpServletRequest request) throws Exception
	{
		String action = request.getParameter("action");

		if (action.equals("removeselect") || action.equals("addselect") || action.equals("toggleselect") || action.equals("selectpage"))
		{
			Set<Long> selectionList = Globals.userSession().selectionMoleculeList;

			if (action.equals("removeselect"))
				selectionList.remove(molecule.id);
			else if (action.equals("addselect") || action.equals("selectpage"))
				selectionList.add(molecule.id);
			else
			{
				// Toggle selection
				if (selectionList.contains(molecule.id))
					selectionList.remove(molecule.id);
				else
					selectionList.add(molecule.id);
			}

		}
		else if ((action.equals("removetag"))  || (action.equals("addtag")))
		{
			Tag tag = 
					assertParam("new-tag")? 
							(Tag) Globals.session().get(Tag.class, getLongParam("new-tag")) :
								(Tag) Globals.session().get(Tag.class, getLongParam("tag"))	;

							if (molecule.tags.contains(tag) && action.equals("removetag"))
								molecule.tags.remove(tag);
							else
								if (!molecule.tags.contains(tag) && action.equals("addtag"))
									molecule.tags.add(tag);

							Globals.session().saveOrUpdate(molecule);
		}
		else 
			throw new Exception("Unknown action - "+action);
	}

	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		WebFilters filters = formFilters(request);
		WebList list = new WebList();

		Criteria mainCriteria = createMainCriteria(filters, Globals.session().createCriteria(Mapping1.class, "mp1")).criteria;
		list.useEntity(Mapping1.class).loadDistinctFromCriteria(mainCriteria, getPageNum(), getPageSize(5));

		for (Object mapObject : list.list) 
		{
			Mapping1 mapping1 = (Mapping1) mapObject;
			//				List<Integer> ids = new ArrayList<Integer>();
			//				for (Mapping2 mapping2 : ((Mapping1) mapping1).mapping2) 
			//					ids.add(mapping2.id);
			//				
			Long count = (Long) Globals.alternateSession().createCriteria(qspr.fragmententities.Mapping1.class)
					.add(Restrictions.eq("id", mapping1.id))
					.add(Restrictions.eq("frag_status", Integer.valueOf(0)))
					.setProjection(Projections.count("id")).uniqueResult();
			if (count > 0)
				mapping1.isFragmented = true;
		}


		return new BrowserModel().setFilters(filters).setObject(list).setTemplate("molecule-browser").getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) 
			throws Exception 
	{		
		int maxAllowed = 50000;
		WebFilters filters = formFilters(request);
		Criteria mainCriteria = createMainCriteria(filters, Globals.session().createCriteria(Mapping1.class, "mp1")).criteria;

		if (!filters.has("id") && filters.in("action", new String[] {"addtag", "removetag", "removeselect"}))
		{
			Set<Long> selectionMoleculeList = Globals.userSession().selectionMoleculeList;

			if (selectionMoleculeList == null || selectionMoleculeList.size() == 0)
				mainCriteria.add(Restrictions.sqlRestriction("1 = 0"));
			else
				mainCriteria.add(Restrictions.in("id", selectionMoleculeList));
		}

		WebList l = new WebList();

		// 1. Select all IDs of matched records
		// 2. Fetch objects for these IDs, by 50 records. 
		//    Restart transaction after each 50-records job
		l.useEntity(Mapping1.class)
		.useProjection(Projections.distinct(Projections.property("id")))
		.loadDistinctFromCriteria(mainCriteria, 1, maxAllowed);

		List<Object> list = l.list;
		if (l.size > maxAllowed)
			throw new UserFriendlyException("You tried to work with " + l.size + " molecules. The maximum allowed number is currently "+maxAllowed+" molecules");

		Iterator<Object> imols = list.iterator();
		while(imols.hasNext())
		{
			int i = 0;
			List<Long> ids = new ArrayList<Long>();
			while (imols.hasNext() && i++ < 50)
				ids.add((Long)imols.next());

			List<Mapping1> molecules = Globals.session().createCriteria(Mapping1.class).add(Restrictions.in("id", ids)).list();
			for (Mapping1 molecule : molecules) 
			{

				logger.info("[EBC] Executing "+request.getParameter("action")+" on record "+molecule.id);
				executeAction((Mapping1)molecule, request);
			}

			Globals.restartAllTransactions(true);
		}

		return new WebModel().getModelAndView();
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) 
			throws Exception
	{
		if(!Globals.userSession().user.isSuperUser())
			throw new UserFriendlyException("Molecules can be accessed only by administrator");
		
		WebModel wm = new WebModel().setList(Globals.getTaginationFilters(Mapping1.class)).setTemplate("molecule-browser");
		return wm.getModelAndView();
	}

	protected void formFilter(String key, String value, WebFilters filters)
	{
		if (key.equals("record-id"))
		{
			ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, getLongParam("record-id"));

			filters.addFilter("id", ep.molecule.mapping1.id.toString(), "");
		} 
		else
			super.formFilter(key, value, filters);
	}
}
