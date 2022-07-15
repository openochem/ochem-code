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

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.entities.Molecule;
import qspr.entities.MoleculeName;
import qspr.entities.ValidatedFact;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.NCBI_Utility;
import qspr.workflow.utils.QSPRConstants;

@Controller
public class NCBISearchController extends BrowserWrapper {

	private static transient final Logger logger = LogManager.getLogger(NCBISearchController.class);

	public NCBISearchController()
	{
		sessionRequired = true;
	}

	@Override
	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws Exception {
		Searcher searcher = (Searcher) request.getSession().getAttribute("searcher");
		return new WebModel(searcher.getWebList()).setTemplate("search-browser").getModelAndView();
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) throws Exception {
		request.getSession().removeAttribute("search-string");
		return new WebModel().setTemplate("search-browser").getModelAndView();
	}

	public ModelAndView searchBrowser(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{  
		// you come here from search in name browser or from the molecule edit search
		Searcher searcher = new Searcher();
		request.getSession().setAttribute("searcher", searcher);

		String searchString = request.getParameter("search-string");
		String emolId = request.getParameter("existmol-id");

		Long existingMolId = null;
		if (emolId != null && ! "".equals(emolId)) {
			existingMolId = Long.parseLong(emolId);
		} 

		searcher.setExistingMoleculeID(existingMolId); // it is important to remember the existing mol id here, otherwise it's lost in the search browser

		// TODO move it up and put the other stuff there also	
		if (null != searchString && ! searchString.equals(""))
			searcher.fetchmolsNCBI(searchString, 1); // changed to search first in pubchem

		return new WebModel().setTemplate("search-browser").getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception {

		String actionParam = request.getParameter("action");

		if ("Check".equals(actionParam))
		{
			logger.info(request.getParameter("Check"));
		}

		if (actionParam.startsWith("Search"))
		{   
			// you come here from search inside the search browser
			Searcher searcher = (Searcher) request.getSession().getAttribute("searcher");
			String searchString = request.getParameter("search-string");

			int searchNumber = 1;
			try {
				searchNumber = Integer.parseInt(request.getParameter("search-number"));
			} catch (NumberFormatException e) {
				// TODO: handle exception
			}

			String dbParam = request.getParameter("database");
			if ("PubChem".equals(dbParam))
			{
				searcher.fetchmolsNCBI(searchString, searchNumber);
			}
			else if ("QSPR".equals(dbParam))
			{
				searcher.fetchQSPRMolecule(searchString, searchNumber);
			}
			System.out.print("");
		}

		return new WebModel().setTemplate("search-browser").getModelAndView();
	}


}

class Searcher {

	private Long existingMoleculeID;

	public Long getExistingMoleculeID() {
		return existingMoleculeID;
	}

	public void setExistingMoleculeID(Long existingMoleculeID) {
		this.existingMoleculeID = existingMoleculeID;
	}

	private List<Molecule> molecules = new ArrayList<Molecule>();

	public void fetchmolsNCBI(String searchString, int resultSize)
	{
		Molecule existingMolecule = (null != existingMoleculeID) ? Repository.molecule.getMolecule(existingMoleculeID) : null;

		try {
			molecules = NCBI_Utility.getMoleculeByName(searchString, resultSize);
			for (Molecule mol : molecules) {
				determineSimilaritySign(mol, existingMolecule);
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	@SuppressWarnings("unchecked")
	public void fetchQSPRMolecule(String searchString, int resultSize)
	{
		//		if (molecules == null)
		//			molecules = new ArrayList<Molecule>();
		molecules.clear();
		Molecule existingMolecule = (null != existingMoleculeID) ? Repository.molecule.getMolecule(existingMoleculeID) : null;

		Criteria criteria = Globals.session().createCriteria(ValidatedFact.class);
		criteria.add(Restrictions.eq("validated", ValidatedFact.VALIDATED))
		.createCriteria("moleculename")
		.add(Restrictions.like("name", searchString));
		criteria.setResultTransformer(Criteria.DISTINCT_ROOT_ENTITY);

		List<ValidatedFact> results = criteria.list();

		for (ValidatedFact vf : results)
		{
			List<String> searchSynonymList = new ArrayList<String>();
			if (vf.mapping != null)
			{
				for (ValidatedFact molVf : vf.mapping.validatedFacts)
				{
					searchSynonymList.add(molVf.moleculename.name);
				}

				List<Molecule> mols = vf.mapping.molecules;
				for (int i=0; i<Math.min(mols.size(), resultSize); i++)
				{
					Molecule m = mols.get(i);

					m.searchSynonymList = searchSynonymList;

					MoleculeName molName = MoleculeName.get(searchString);
					if (molName != null)
						m.searchedBy = molName;

					String smiles = "";
					try
					{
						smiles = Various.molecule.convertToSmilesOrSmart(m.getData(),QSPRConstants.SMILESNOAROM);
					} catch (IOException e)
					{
						e.printStackTrace();
					}

					m.smile = smiles;

					determineSimilaritySign(m, existingMolecule);
					molecules.add(m);

				}
			}			
		} 
	}

	private void determineSimilaritySign(Molecule foundMol, Molecule existMol) 
	{ 
		if (null == existMol)
			foundMol.searchCompSign = "!=";
		else {

			foundMol.searchExistingMolId = existingMoleculeID;



			if (foundMol.mapping1.id.intValue() == existMol.mapping1.id.intValue()) {
				if (foundMol.mapping2.id.longValue() == existMol.mapping2.id.longValue()) {
					foundMol.searchCompSign = "=";
				} else {
					//					String almostEqual = (new Character((char)248)).toString();
					foundMol.searchCompSign = "~"; //almostEqual;
				}
			}
			else { 
				foundMol.searchCompSign = "!=";
			}
		}
	}

	public WebList getWebList() {
		WebList webList = new WebList();
		return webList.loadFromList(molecules); 
	}

}
