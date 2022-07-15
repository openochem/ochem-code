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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.business.WebFilters;
import qspr.dao.Repository;
import qspr.entities.Mapping2;
import qspr.entities.Molecule;
import qspr.frontend.MoleculeProfileData;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.util.MoleculePeer;
import qspr.util.NCBI_Utility;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.descriptorcache.CacheEntry;
import com.eadmet.descriptorcache.DescriptorConfigEntry;
import com.eadmet.descriptorcache.DescriptorsRepository;
import com.eadmet.descriptorcache.DescriptorsRepositoryFactory;
import com.eadmet.exceptions.UserFriendlyException;

@Controller
public class MoleculeController extends BrowserWrapper 
{
	public MoleculeController()
	{
		//		Commented out because in e.g. OchemTrade context the getFromText function has to be called with no login requirements. 
		//		Alternate solution - separate that function into separate controller (resort to it if current solution causes obvious problems) Nos 21.08.2014
		//		sessionRequired = true;

	}

	/**
	 * Display the molecule profile page
	 */
	public ModelAndView profile(HttpServletRequest request, HttpServletResponse response) throws IOException
	{
		if (Globals.userSession() == null || Globals.userSession().user == null)
			return redirect("user/pleaseregister.do");
		
		if(!Globals.isValidatedUser())
			throw new UserFriendlyException(QSPRConstants.RESTRICT);
		
		Mapping2 mp2;
		if (assertParam("depiction"))
			mp2 = Repository.molecule.getMolecule(getLongParam("depiction")).mapping2;
		else
			mp2 = Repository.molecule.getMapping2(getIntParam("id"));
		MoleculeProfileData mpData = new MoleculeProfileData(mp2);

		return new WebModel(mpData).setTemplate("molecule-profile").getModelAndView();
	}

	/**
	 * Get the descriptors storage entries for a particular molecule
	 */
	public ModelAndView getDescriptors(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		int mp2 = getIntParam("mp2");

		DescriptorsRepository dRepository = DescriptorsRepositoryFactory.getReattemptingRepository();

		List<CacheEntry> cacheEntries = new ArrayList<CacheEntry>();

		for (String user : new String[]{Globals.userSession().user.login, null})
		{
			List<DescriptorConfigEntry> dConfigs = dRepository.getConfigurations(user);
			for (DescriptorConfigEntry dConfig : dConfigs)
			{
				dConfig.setUser(user);
				List<CacheEntry> entries = dRepository.getDescriptors(new Integer[]{mp2}, dConfig);
				if (entries.get(0) != null)
					cacheEntries.add(entries.get(0));
			}
		}

		return new WebModel().setList(cacheEntries).getModelAndView();
	}

	public ModelAndView action(HttpServletRequest request, HttpServletResponse response) throws Exception{

		if (request.getParameter("action").equals("Get"))
		{
			return internalGet(request);
		} else
			if (request.getParameter("action").equals("Submit") 
					|| request.getParameter("action").equals("Update") 
					|| request.getParameter("action").equals("Smiles")
					|| request.getParameter("action").equals("Searcher"))
			{
				return internalSubmit(request);
			} 
		return new WebModel().getModelAndView();
	}

	@SuppressWarnings("unchecked")
	private ModelAndView internalGet(HttpServletRequest request) throws Exception
	{
		Molecule molecule;
		String str = request.getParameter("mol-id");
		if (str.matches("[0-9]+"))
			molecule = Repository.molecule.getMolecule(Long.valueOf(str));
		else
			molecule = MoleculePeer.fetchFromString(str, null);
		List<Object> similar = null;
		if (molecule != null)
		{
			similar = Globals.session().createCriteria(Molecule.class).add(Restrictions.eq("mapping1", molecule.mapping1)).list();
		}
		return new WebModel(molecule).setList(similar).getModelAndView();
	}

	public ModelAndView getFromText(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Molecule mol = MoleculePeer.fetchFromString(getParam("molecule-str"), null);
		return new WebModel(mol).getModelAndView();
	}

	@SuppressWarnings("unchecked")
	private ModelAndView internalSubmit(HttpServletRequest request) throws Exception
	{
		Molecule molecule = null;

		File f = null;
		try
		{
			f = Globals.getUploadedFile();
		} catch (Exception e)
		{
			//e.printStackTrace();
		}

		if (f != null)
		{
			molecule = MoleculePeer.getMolecule(f);
		}
		else
		{
			if ("Searcher".equals(request.getParameter("action"))){
				Molecule pubChemMol = NCBI_Utility.getMoleculeByName(request.getParameter("Search"));
				if (pubChemMol != null){
					request.getSession().setAttribute("search-mol", pubChemMol);
					molecule = pubChemMol;
				}
			} 
			else if(assertParam("mol-id"))
				molecule = Repository.molecule.getMolecule(Long.valueOf(request.getParameter("mol-id")));
			else if(assertParam("data"))
				molecule = MoleculePeer.fetchFromString(request.getParameter("data"), null);
			//molecule = MoleculePeer.getMolecule(request.getParameter("data"));

		}

		List<Object> similar = null;
		if (molecule != null)
		{
			similar = Globals.session().createCriteria(Molecule.class).add(Restrictions.eq("mapping1", molecule.mapping1)).list();
		}
		return new WebModel(molecule).setList(similar).getModelAndView();
	}

	public ModelAndView show(HttpServletRequest request, HttpServletResponse response) 
			throws Exception
	{
		WebModel wm = new WebModel().setTemplate("molecule");
		return wm.getModelAndView();
	}

	@SuppressWarnings("unchecked")
	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		Molecule molecule = null;
		WebFilters filters = formFilters(request);

		//molecule editor
		if (filters.has("mol-id") && filters.getLong("mol-id") > 0){
			molecule = Repository.molecule.getMolecule(Long.valueOf(filters.get("mol-id")));

			if (filters.has("browser") && filters.get("browser").equals("depictions"))
			{
				List<Object> similar = Globals.session()
						.createCriteria(Molecule.class)
						.add(Restrictions.eq("mapping1", molecule.mapping1))
						.list();
				return new WebModel(new WebList().loadFromList(similar)).getModelAndView(); 

			} else 
			{
				if (!assertParam("format"))
					molecule = Repository.molecule.getMolecule(Long.valueOf(request.getParameter("id")));

				List<Object> molecules = Globals.session().createCriteria(Molecule.class).add(Restrictions.eq("mapping1", molecule.mapping1)).list();
				return new WebModel(new WebList().loadFromList(molecules)).getModelAndView();
			}
		}
		return null;
	}

	public ModelAndView edit(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		WebModel wm = new WebModel().setTemplate("molecule");
		return wm.getModelAndView();
	}

}
