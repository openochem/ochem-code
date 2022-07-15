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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.modelling.applier.PropertyPrediction;

import com.eadmet.utils.DoubleMap;
import com.eadmet.utils.NumericalValueStandardizer;

/**
 * A subset of all MMPs. For example, MMPs for a selected property or MMPs for molecules form a particular basket.
 * The MMPs are grouped by transformations. 
 * 
 * Since the subsets are cached in memory, only IDs of the MMPs are stored.
 * 
 * @author midnighter
 *
 */
public class MMPSubset 
{
	public String id;
	public long pairsCount = 0;

	private LinkedHashMap<Long, List<ThinPair>> pairsByTransformations = new LinkedHashMap<Long, List<ThinPair>>();
	public Map<Integer, PropertyPrediction[]> predictions = new HashMap<Integer, PropertyPrediction[]>();

	/**
	 * Experimental values for the molecules within this subset (if available)
	 */
	public DoubleMap<Long, Integer, FuzzyValue> molValues = new DoubleMap<Long, Integer, FuzzyValue>();

	public List<ThinPair> pairs = new ArrayList<ThinPair>();

	public void addPair(Long transformationId, Long mmpId, int mol1, int mol2) {
		addPair(new ThinPair(mmpId, mol1, mol2, transformationId));
	}

	public void addPair(ThinPair pair) {
		List<ThinPair> trPairs = pairsByTransformations.get(pair.transformationId);
		if (trPairs == null)
			pairsByTransformations.put(pair.transformationId, trPairs = new ArrayList<ThinPair>());

		trPairs.add(pair);
		pairs.add(pair);
		pairsCount++;
	}

	/**
	 * A memory- and time- efficient way to load large number of molecular matched pairs by their identifiers
	 */
	public List<ThinPair> getThinPairs(Long transformationId) {

		if (transformationId != null)
			return pairsByTransformations.get(transformationId);
		else
			return pairs;
	}

	/**
	 * A memory- and time- efficient way to load large number of molecular matched pairs by their identifiers
	 */
	public List<MMPair> getPairs(Long transformationId) {
		List<MMPair> pairs = new ArrayList<MMPair>();

		for (ThinPair pair : getThinPairs(transformationId))
			pairs.add(MMPair.thin(pair.id, pair.mol1Id, pair.mol2Id));

		return pairs;
	}


	public void fillPairsWithData(List<MMPair> pairs) {
		if (molValues.isEmpty())
			return;
		Property property = Property.getById(molValues.keySet().iterator().next());
		Map<Integer, FuzzyValue> vals = molValues.get(property.id);
		for (MMPair pair : pairs)
			if (vals.containsKey(pair.molecule1.id) && vals.containsKey(pair.molecule2.id))
				if (property.isQualitative())
				{
					PropertyOption option1 = PropertyOption.getById(vals.get(pair.molecule1.id).getLongValue());
					PropertyOption option2 = PropertyOption.getById(vals.get(pair.molecule2.id).getLongValue());
					pair.values = new String[]{option1.name, option2.name};
				}
				else
					pair.values = new String[]{NumericalValueStandardizer.getSignificantDigits(vals.get(pair.molecule1.id).value), NumericalValueStandardizer.getSignificantDigits(vals.get(pair.molecule2.id).value)};
	}

	public Set<Long> getTransformationIds() {
		return pairsByTransformations.keySet();
	}

	public List<Long> getPairIds() {
		List<Long> ids = new ArrayList<Long>();
		for (ThinPair pair : pairs)
			ids.add(pair.id);

		return ids;
	}

	@SuppressWarnings("unchecked")
	public Set<Integer> getMoleculeIds() {
		Set<Integer> uniqueMols = new HashSet<Integer>();
		List<Long> pairIds = getPairIds();

		if (pairIds.size() == 0)
			return uniqueMols;

		List<Integer> mols = Globals.session().createCriteria(MMPair.class).add(Restrictions.in("id", pairIds))
				.setProjection(Projections.property("molecule1")).list();
		uniqueMols.addAll(mols);

		mols = Globals.session().createCriteria(MMPair.class).add(Restrictions.in("id", pairIds))
				.setProjection(Projections.property("molecule2")).list();
		uniqueMols.addAll(mols);

		return uniqueMols;
	}

	public MMPSubset clone() {
		MMPSubset clone = new MMPSubset();

		for (ThinPair pair : pairs)
			clone.addPair(pair);

		clone.molValues = molValues.clone();

		return clone;
	}

}
