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

package qspr.business;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.HibernateException;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.entities.Mapping2Filter;

import com.eadmet.utils.OCHEMUtils;


/**
 * The filters related to substructure/similarity search
 * 
 * @author midnighter
 *
 */
public class SubstructureSearchFilter
{
	/**
	 * The structure (SMILES or SDF) of the query molecule
	 */
	public String structure;

	/**
	 * The SMARTS pattern for substructure query
	 */
	public String smarts;


	/**
	 * The saved set (created for long queries)
	 */
	public Integer setid;

	/**
	 * Maximum number of retrieved molecules; can be modified by the user
	 */
	public Integer maximumHits = 1000;

	/**
	 * The search type (similarity or substructure)
	 */
	public String type;

	/**
	 * The threshold of the similarity search (in percents)
	 */
	public double similarityThreshold;

	private String smiles;

	/**
	 * The cache of substructure filters. Maps md5(structure + search type) to the filter ID
	 */
	private static Map<String, Integer> substructureFiltersCache = new HashMap<String, Integer>();

	public SubstructureSearchFilter()
	{

	}

	/**
	 * Is it a substructure search? (alternatively, it might me a similarity search)
	 * @return
	 */
	public boolean isSubstructureSearch()
	{
		return "substructure".equals(type);
	}

	/**
	 * Get the SMILES representation of the queried structure
	 */
	private String getSmiles() throws IOException
	{
		if (smiles != null)
			return smiles;

		return smiles = Various.molecule.convertToSmilesOrSmart(structure);
	}


	/**
	 * Get a Mapping2 filter ID corresponding to the matched molecules
	 * Currently not used since global search is disabled
	 */
	public int getFilterID() throws HibernateException, IOException
	{
		String md5 = OCHEMUtils.getMD5(getSmiles() + type + similarityThreshold);
		if (substructureFiltersCache.containsKey(md5))
			return substructureFiltersCache.get(md5);
		else
		{
			int filterID = Mapping2Filter.generateFilterID();
			logger.info("Generating a new substructure search filter for " + getSmiles() +", ID " + filterID);
			if (isSubstructureSearch())
				Globals.session().createSQLQuery("insert into Mapping2Filter(filter_id, mapping2_id) select " + filterID + ", mapping2_id from StructureQuery where match_substructure(screen, molecule, '" + getSmiles() + "')").executeUpdate();
			else
				Globals.session().createSQLQuery("insert into Mapping2Filter(filter_id, mapping2_id) select " + filterID + ", mapping2_id from StructureQuery where similarity('" + getSmiles() + "', simscreen) > " + similarityThreshold).executeUpdate();
			substructureFiltersCache.put(md5, filterID);

			return filterID;
		}
	}

	/**
	 * Construct filters from user's input via HTTP
	 * @param request
	 */
	public SubstructureSearchFilter(WebFilters filters)
	{
		if (filters.has("xemistry-smiles"))
		{
			structure = filters.get("xemistry-smiles").trim();
			if (structure.matches("M[0-9]+"))
				structure = Repository.molecule.getMapping2(Integer.valueOf(structure.substring(1))).getMolecule().getData();
			else if (structure.matches("DP[0-9]+"))
				structure = Repository.molecule.getMolecule(Integer.valueOf(structure.substring(2))).getData();
			type = filters.get("xemistry-sstype");
			similarityThreshold = Double.valueOf(filters.get("xemistry-similarity-cutoff"));
		} else
			if (filters.has("xemistry-smarts"))
				smarts = filters.get("xemistry-smarts").trim();

		if (filters.has("saved-set-id"))
			setid = Integer.parseInt(filters.get("saved-set-id"));

		if (filters.has("maximum-hits"))
			maximumHits = Integer.parseInt(filters.get("maximum-hits"));

	}

	private static final Logger logger = LogManager.getLogger(SubstructureSearchFilter.class);
}
