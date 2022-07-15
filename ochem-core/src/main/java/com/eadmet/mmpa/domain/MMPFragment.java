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

package com.eadmet.mmpa.domain;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;

import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.dao.Various;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.MemoryUtils;

@Entity
public class MMPFragment
{
	@Id
	@Column(name = "fragment_id")
	@GeneratedValue
	public Long id;

	@Column
	public String inchi;

	@Column
	public String smiles;

	/**
	 * Number of atoms
	 */
	@Column
	public byte size;

	/**
	 * Is this fragment a simple unbranched chain of carbons?
	 */
	@Column(name = "carbon_chain")
	public Boolean carbonChain;

	private static Map<String, MMPFragment> cache = new HashMap<String, MMPFragment>();
	private static final Logger logger = LogManager.getLogger(MMPFragment.class);

	static public void clearCache() {
		cache.clear();
	}

	public static MMPFragment getFragment(Long id) {
		return (MMPFragment) Globals.session().get(MMPFragment.class, id);
	}

	public static MMPFragment getFragment(String inchi, String smiles) {

		if (cache.isEmpty()) {
			logger.info("Filling the MMPFragment cache");
			logger.info(MemoryUtils.memorySummary());
			@SuppressWarnings("unchecked")
			List<MMPFragment> frags = Globals.session().createCriteria(MMPFragment.class).list();
			for (MMPFragment frag : frags)
				cache.put(frag.inchi, frag);
			logger.info("MMPFragment cache filled with " + cache.size() + " elements");
			logger.info(MemoryUtils.memorySummary());
		}

		MMPFragment frag = cache.get(inchi);

		if (frag == null) {
			frag = new MMPFragment();
			frag.inchi = inchi;
			frag.smiles = smiles;
			frag.smilesUpdated();
			Globals.session().save(frag);
			cache.put(inchi, frag);
		}

		return frag;
	}

	public void smilesUpdated() {
		carbonChain = smiles.replaceAll("\\[Al\\]", "").matches("C+");

		try {
			int atomCount = Various.molecule.getAtomCount(Various.molecule.convertToFormat(smiles, QSPRConstants.SDFNOH));
			int alCount = StringUtils.countMatches(smiles, "Al");
			size = Byte.valueOf("" + (atomCount - alCount));
		} catch (NumberFormatException e) {
			size = 0;
		} catch (IOException e) {
			size = 0;
			logger.warn("A problem with fragment " + smiles);
			size = 0;
			e.printStackTrace();
		}
		
	}

}
