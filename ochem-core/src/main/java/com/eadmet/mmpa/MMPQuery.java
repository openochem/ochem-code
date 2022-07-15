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

package com.eadmet.mmpa;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.apache.commons.lang.StringUtils;

import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.Property;
import qspr.modelling.applier.ModelApplier;

public class MMPQuery
{
	public Long basketId;
	public Long propertyId;
	public Long tagId;
	public Long transformationId;
	
	public Model model;
	
	public boolean significantTransformationsOnly;
	public boolean affectedPairsOnly;
	
	/**
	 * null = no additional filtering, 1 = only pairs with increasing property, -1 = only pairs with decreasing property
	 * works only when requestData (below) is set, and transformationId is not null, and propertyId is not null
	 */
	
	public Integer propertyChangeDirection;
	
	/**
	 * Also request the experimental values for the molecules in pairs
	 */
	public boolean requestData;
	
	/**
	 * Minimum Tanimoto similarity netween selected pairs (0-100)
	 */
	public Short similarity;
	
	/**
	 * Molecular filter ID
	 */
	public Integer filterId;
	
	public Collection<Integer> molIDs;
	
	public ModelApplier applier;
	
	public boolean insideAD;
	
	/**
	 * Maximum number of allowed pairs. If more available, try limiting by similarity.
	 * Sometimes we can't afford to query/display too many pairs (for example, in UI). 
	 * This is what this parameter is for.
	 */
	public Long maxPairs;
	
	public MMPQuery setBasket(Long basketId)
	{
		this.basketId = basketId;
		return this;
	}
	
	public MMPQuery setProperty(Long propertyId)
	{
		this.propertyId = propertyId;
		return this;
	}
	
	public MMPQuery setSimilarity(Short similarity)
	{
		this.similarity = similarity;
		return this;
	}
	
	public String toString() {
		List<String> res = new ArrayList<String>();
		if (basketId != null) 
			res.add("Basket " + Basket.getById(basketId).name);
		
		if (propertyId != null) 
			res.add("Property " + Property.getById(propertyId).getName());
		
		if (transformationId != null)
			res.add("Transformation " + transformationId);
		
		if (filterId != null)
			res.add("MolFilter " + filterId);
		
		if (similarity != null)
			res.add("Similarity " + similarity);
		
		return StringUtils.join(res, ", ");
	}
	
	public boolean hasFilters() {
		return basketId != null || propertyId != null || transformationId != null || filterId != null || tagId != null;
	}
	
	public MMPQuery() {
		
	}
	
	public MMPQuery(Basket b) {
		this.basketId = b.id;
	}
	
	public MMPQuery(ModelApplier applier, boolean insideAD) {
		this.applier = applier;
		this.insideAD = insideAD;
	}

}
