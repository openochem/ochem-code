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
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Molecule;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.CASRN;
import qspr.util.MoleculePeer;
import qspr.util.NCBI_Utility;

@Controller
public class MoleculeSearchController extends ControllerWrapper 
{
	private static transient final Logger logger = LogManager.getLogger(MoleculeSearchController.class);
	
	public MoleculeSearchController()
	{
		sessionRequired = true;
	}
	
	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		String moleculeString = getParam("name");
		List<Molecule> results = new ArrayList<Molecule>();

		if (moleculeString.matches("M[0-9]+")) // MoleculeID
		try
		{
			results.add(MoleculePeer.getByID2(moleculeString));
		} catch (Exception e)
		{
			logger.debug(e);
		}
		
		if (moleculeString.matches("BM[0-9]+")) // MoleculeID
		try
		{
			results.add(MoleculePeer.getByID1(moleculeString));
		} catch (Exception e)
		{
			logger.debug(e);
		}
			
		if (moleculeString.equals("empty"))
			results.add(Repository.molecule.getEmptyMolecule());
			
		try
		{
			results.add(MoleculePeer.getByStructure(moleculeString));
		} catch (Exception e)
		{
			//Not a valid molecule structure, continue
		}
			
		String name = CASRN.checkCasrnSyntax(moleculeString);
		if (name == null) //Not CASRN
			name = moleculeString;
		
		results.addAll(MoleculePeer.getByValidatedFact(moleculeString));	
		
		if (results.size() > 0) //Only do PubChem stuff if necessary
			return new WebModel(new WebList().loadFromList(results)).getModelAndView();
			
		Globals.restartAllTransactions(true);
		results = NCBI_Utility.getMoleculeByName(name, 10);
		Globals.restartAllTransactions(true);
		
		return new WebModel(new WebList().loadFromList(results)).getModelAndView();
	}
}
