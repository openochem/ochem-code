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

import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.Criteria;
import org.hibernate.Hibernate;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.MarshallingOption;
import qspr.dao.Repository;
import qspr.entities.Alert;
import qspr.entities.Article;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Molecule;
import qspr.entities.MoleculeName;
import qspr.entities.Predicate;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.PropertyValue;
import qspr.entities.Unit;
import qspr.exception.DublicateRecordException;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.AccessChecker;
import qspr.util.MoleculePeer;

import com.eadmet.business.MoleculeNamesService;
import com.eadmet.business.PaginationFilter;
import com.eadmet.exceptions.UserFriendlyException;

// This class was previously a part of EPBrowserController
// And appeared as a separated class as a result of refactoring
// Midnighter

@Controller
public class EPRecordController extends BrowserWrapper
{
	MoleculeNamesService mnService = new MoleculeNamesService();
	public EPRecordController()
	{
		sessionRequired = true;
	}

	@Override
	public ModelAndView list(HttpServletRequest req, HttpServletResponse res)
			throws Exception 
	{
		return null;
	}

	@SuppressWarnings("unchecked")
	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		ExperimentalProperty ep = new ExperimentalProperty();
		Molecule m = null;

		if(assertParam("id") && getLongParam("id") > 0)
		{
			ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, getLongParam("id"));
			AccessChecker.requestModificationPermission(ep);
		}

		if (assertParam("hide"))
			ep.rights = Globals.RIGHTS_NONE;
		else
			ep.rights = Globals.RIGHTS_FREELY_AVAILABLE;
		//update all experimental records of same molecule within a article if condition is true
		//Anil
		if (assertParam("n-molecule"))
			m = Repository.molecule.getMolecule(Long.valueOf(request.getParameter("n-molecule")));

		if((assertParam("modify-similar-ep")))
		{
			Criteria criteria =  Globals.session().createCriteria(ExperimentalProperty.class)
					.createAlias("molecule", "mol")
					.add(Restrictions.eq("article", ep.article))
					.add(Restrictions.isNull("deleted"));
			if(ep.moleculenames.size() > 0){
				criteria.createCriteria("moleculenames")
				.add(Restrictions.like("name", "%"+ep.moleculenames.get(0)+"%"));
			}
			List<ExperimentalProperty> similar_ep = criteria.add(Restrictions.eq("mol.mapping1", ep.molecule.mapping1)).list();

			for (ExperimentalProperty experimentalProperty : similar_ep) {
				if(m != null)
				{
					experimentalProperty.molecule = m;
					Globals.session().saveOrUpdate(experimentalProperty);
				}
			}
		}

		if (assertParam("n-property"))
			ep.property = (Property) Globals.session().get(Property.class, Long.valueOf(request.getParameter("n-property")));
		if (assertParam("n-value"))
			ep.setValueWithPredicate(request.getParameter("n-value"));
		if (assertParam("n-predicate"))
			ep.predicate = (Predicate) Globals.session().get(Predicate.class, Long.valueOf(request.getParameter("n-predicate"))); 
		if(assertParam("n-second-value"))
		{
			ep.secondValue = Double.valueOf(request.getParameter("n-second-value"));
			if (ep.secondValue < 0 && "+-".equals(ep.predicate.shortName))
				throw new UserFriendlyException("Accuracy should pe a positive number. You provided "+ep.secondValue);
			if (ep.secondValue < ep.value && "-".equals(ep.predicate.shortName))
				throw new UserFriendlyException("Left bound of interval should be less then right one.");
		}
		else
			ep.secondValue = null;

		if (assertParam("n-unit"))
			ep.unit = (Unit) Globals.session().get(Unit.class, Long.valueOf(request.getParameter("n-unit")));
		if (assertParam("n-molecule"))
		{   
			//Inchi check irrelevant here, even if inchi is the same, we should set the new molecule
			ep.molecule = m;
		} else
		{
			//Provide empty molecule;
			ep.molecule = MoleculePeer.getMolecule("");
		}

		ep.artPageNum = getIntParam("n-page");
		ep.artLineNum = getIntParam("n-line");
		ep.artTableNum = request.getParameter("n-table");
		ep.artMolId = request.getParameter("n-art-mol-id");
		ep.other = request.getParameter("n-comment");

		if (assertParam("n-article"))
			ep.article = (Article) Globals.session().get(Article.class, getLongParam("n-article"));

		if (ep.article == null)
			throw new UserFriendlyException("You must specify an article to save the record");

		if (assertParam("n-conditions"))
			ep.conditions = (ConditionSet) Globals.session().get(ConditionSet.class, getLongParam("n-conditions"));
		if (assertParam("n-option"))
			ep.setOption((PropertyOption) Globals.session().get(PropertyOption.class, getLongParam("n-option")));

		// ep.status is used to put the status in evidence for error in record(0) or not validated(1)
		if (assertParam("evidence"))
		{
			switch (getIntParam("evidence"))
			{
			case 0: // no evidence
			{
				ep.connectedProperty = null; 
				ep.ep_status = null;
				break;
			}
			case 1: // measured in this article 
			{
				ep.connectedProperty = ep; 
				ep.ep_status = null;
				break;
			}
			case 2: // measured in another article
			{
				ep.connectedProperty = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, getLongParam("n-reference"));
				ep.ep_status = null;
				break;
			}
			case 3: // measured in this article (to be verified)
			{
				ep.connectedProperty = ep; 
				ep.ep_status = ExperimentalProperty.STATUS_TOVERIFY;
				break;
			}
			case 4: // error in the record (see discussion)
				ep.ep_status = ExperimentalProperty.STATUS_ERROR;
				break;
			case 5: // invalid record
				ep.ep_status = ExperimentalProperty.STATUS_INVALID;
				break;
			}
		}

		if (ep.property.isQualitative() && ep.option == null)
			throw new UserFriendlyException("You did not select value for property "+ep.property.getName());

		// Save conditions. Now record and its conditions are saved during the same request / Midnighter
		String[] conditionIds = request.getParameterValues("cond-id");
		String[] conditionValues = request.getParameterValues("cond-value");
		String[] conditionTextValues = request.getParameterValues("cond-text-value");
		String[] conditionOptions = request.getParameterValues("cond-option");
		String[] conditionUnits = request.getParameterValues("cond-unit");
		String[] conditionPredicates = request.getParameterValues("cond-pred");
		String[] conditionSecondValues = request.getParameterValues("cond-second-value");
		ConditionSet conditionSet = new ConditionSet();
		if (conditionIds != null)
			for (int i = 0; i < conditionIds.length; i++)
			{
				PropertyValue cv = new PropertyValue();
				cv.property = (Property) Globals.session().get(Property.class, Long.valueOf(conditionIds[i]));
				if (cv.property.isQualitative())
				{
					if (!"".equals(conditionOptions[i]))
						cv.option = (PropertyOption) Globals.session().get(PropertyOption.class, Long.valueOf(conditionOptions[i]));
					assert cv.option != null;
				}
				else if (cv.property.isNumeric())
				{
					Double val = Double.valueOf(conditionValues[i]);
					cv.unit = (Unit) Globals.session().get(Unit.class, Long.valueOf(conditionUnits[i]));
					cv.value = val;
					if (conditionPredicates[i] != null && !conditionPredicates[i].equals(""))
						cv.predicate = (Predicate) Globals.session().get(Predicate.class, Long.valueOf(conditionPredicates[i])); 
					if (conditionSecondValues[i] != null && !conditionSecondValues[i].equals(""))
					{
						cv.secondValue = Double.valueOf(conditionSecondValues[i]);
						if (cv.secondValue < 0 && "+-".equals(cv.predicate.shortName))
							throw new UserFriendlyException("Accuracy should pe a positive number. You provided "+cv.secondValue+" for condition "+cv.property.getName());
						if (cv.secondValue < cv.value && "-".equals(cv.predicate.shortName))
							throw new UserFriendlyException("Left bound of interval should be less then right one for condition "+cv.property.getName()+".");
					}
					else
						cv.secondValue = null;

					assert cv.unit != null;
				}
				else
					cv.textualValue = conditionTextValues[i];

				conditionSet.values.add(cv);
			}
		ep.conditions = conditionSet.get();

		// Ownership and rights		
		if (Globals.userSession().user != null)
		{
			if (ep.id == null || ep.id < 0)
				ep.introducer = Globals.userSession().user;

			ep.owner = Globals.userSession().user;
			if (assertParam("rights"))
			{
				ep.rights = Integer.valueOf(request.getParameter("rights"));
				Globals.userSession().user.defaultRights = ep.rights;
			}
		}

		ep.updateHash();
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
				throw new DublicateRecordException(ep.duplicate);

		Globals.session().saveOrUpdate(ep);

		// end of edit
		return new WebModel(ep).getModelAndView();
	}

	@SuppressWarnings("unchecked")
	public ModelAndView edit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		if (Globals.userSession().user == null)
			return redirect("user/pleaseregister.do");

		Globals.setMarshallingOption(MarshallingOption.PROPERTY_OPTIONS);
		ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, getLongParam("id"));
		if (ep != null)
		{
			Criteria criteria = Globals.session().createCriteria(ExperimentalProperty.class)
					.createAlias("molecule", "mol")
					.add(Restrictions.eq("article", ep.article))
					.add(Restrictions.isNull("deleted"));
			if(ep.moleculenames.size() > 0)
			{
				criteria.createCriteria("moleculenames")
				.add(Restrictions.like("name", "%"+ep.moleculenames.get(0)+"%"));
			}
			criteria
			.add(Restrictions.eq("mol.mapping1", ep.molecule.mapping1)).setProjection(Projections.rowCount());

			ep.count = ((Long)criteria.list().get(0)).intValue();
		}
		else{
			ep = new ExperimentalProperty();
			ep.id = Long.valueOf(-1);
			if (assertParam("clone"))
			{
				ExperimentalProperty source = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, getLongParam("clone"));
				//Hibernate.initialize(source.getCheckedNames());
				Hibernate.initialize(source.moleculenames);
				Hibernate.initialize(source.colorednames);
				Globals.session().evict(source);
				source.id = Long.valueOf(-1);
				source.introducer = Globals.userSession().user;
				source.connectedProperty = null;
				//source.basketEntries = null;
				ep = source;
			}
		}
		List<Predicate> predicates = Globals.session().createCriteria(Predicate.class).list();
		return new WebModel(ep).setList(predicates).setTemplate("record-edit").getModelAndView();
	}

	public ModelAndView nameactions(HttpServletRequest request, HttpServletResponse res) throws Exception
	{
		if ("checknames".equals(getParam("action")))
		{
			Long ep_id = getLongParam("id");
			if (ep_id != null &&  ep_id != -1)
			{
				ExperimentalProperty ep = (ExperimentalProperty)Globals.session().get(ExperimentalProperty.class, ep_id);
				mnService.checkNames(ep);
			}
			return new WebModel().getModelAndView();
		}
		else if ("validate".equals(getParam("action"))) //name validated by user
		{
			if(Globals.userSession().user == null)
				throw new UserFriendlyException("Guest users are not allowed to validate molecule names");

			Long mn_id = getLongParam("id");
			Long mol_id = getLongParam("mol-id");
			if (mn_id != null && mol_id != null)
			{
				MoleculeName mn = (MoleculeName) Globals.session().get(MoleculeName.class, Long.valueOf(request.getParameter("id")));
				Molecule m = Repository.molecule.getMolecule(Long.valueOf(request.getParameter("mol-id")));
				mnService.validate(mn, m);
			}

			return new WebModel(new Alert("ok")).getModelAndView();
		}
		else if(request.getParameter("action").equals("invalidate") && (assertParam("id"))) //name invalidated by user
		{
			if(Globals.userSession().user == null)
				throw new UserFriendlyException("Guest users are not allowed to invalidate molecule names");
			Long mn_id = getLongParam("id");
			if (mn_id != null)
			{
				MoleculeName mn = (MoleculeName) Globals.session().get(MoleculeName.class, Long.valueOf(request.getParameter("id")));
				mnService.invalidate(mn);
			}
			return new WebModel(new Alert("ok")).getModelAndView();
		}
		else  //Save names for the ExperimentalProperties
		{  
			Long ep_id = getLongParam("id");
			if (ep_id != null && ep_id != -1)
			{
				ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, getLongParam("id"));
				String[] names = request.getParameterValues("moleculename");
				mnService.setNamesForEP(ep, names);
			}
			return new WebModel(new Alert("ok")).getModelAndView();
		}
	}

	public ModelAndView listnames(HttpServletRequest request, HttpServletResponse response)
	{
		WebList wl = new WebList();
		Long ep_id = getLongParam("id");
		if (ep_id != null && ep_id != -1)
		{
			ExperimentalProperty ep = (ExperimentalProperty)Globals.session().get(ExperimentalProperty.class, ep_id);
			Molecule mol = null;
			Long mol_id = getLongParam("mol-id");

			if (mol_id != null)
				mol = Repository.molecule.getMolecule(mol_id);

			wl.loadFromList(mnService.listMoleculeNames(ep, mol));
			Globals.session().evict(ep);
		}
		return new WebModel(wl).getModelAndView();
	}

	public ModelAndView getsynonyms(HttpServletRequest req, HttpServletResponse res)
	{	
		WebList wl = new WebList();

		Long molId = getLongParam("mol-id");

		if (molId == null)
			return new WebModel(wl).getModelAndView();

		if (molId > 0)
		{
			Molecule mol = Repository.molecule.getMolecule(molId);

			Long epId = getLongParam("id");
			ExperimentalProperty ep = null;

			if (epId != null)
				ep =  (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, epId);

			PaginationFilter pager = new PaginationFilter(getPageNum(), getPageSize(10));
			List<MoleculeName> names = mnService.listSynonyms(mol, ep, pager);
			wl.loadFromPagedList(names, pager.pageNum, pager.pageSize, pager.totalSize);
		} else 
		{
			Molecule mol = (Molecule) req.getSession().getAttribute("search-mol");
			List<MoleculeName> names = mnService.listSynonymsFromSearch(mol);
			wl.loadFromList(names);
		}
		return new WebModel(wl).getModelAndView();
	}

}
