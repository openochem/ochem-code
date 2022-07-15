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

import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang.StringUtils;

import com.eadmet.mmpa.domain.MMPSubset;

import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.Property;
import qspr.modelling.applier.ModelApplier;

public class MMPStatsRequest
{
	public Model model;
	
	/**
	 * When calculating statistics based on a model, use the predicted values
	 */
	public boolean usePredictedValues;
	
	public Basket basket;
	
	public Property property;
	
	public ModelApplier applier;
	
	/**
	 * The largest acceptable pValue to filter "significant" transformations.
	 * "1" stands for no pValue filtering
	 */
	public double pValue = 1;
	
	public boolean insideAD;
	
	/**
	 * This is a fancy one. 
	 * If set, the statistics will be bootstrapped taking into account the estimated prediction accuracy.
	 */
	public boolean bootstrapStatistics;
	
	public int minPairs = 4;
	
	public Integer maxAtoms;
	
	public Double meanStdFactor;
	
	public MMPSubset subset;
		
	public String transformationPattern;
	
	public boolean publicData;
	public boolean primaryRecords;
	
	public Short similarity;
	
	public Long maxPairs;
	
	public String getMinorKey() {
		StringWriter sw = new StringWriter();
		if (model != null)
			sw.append("-m-" + model.id + "-" + usePredictedValues);
		if (basket != null)
			sw.append("" + basket);
		if (similarity != null)
			sw.append("-sim-" + similarity);
		if (property != null)
			sw.append("-p-" + property.id);
		if (applier != null)
			sw.append("-a-" + applier);
		if (bootstrapStatistics)
			sw.append("-bootstrap");
		if (subset != null)
			sw.append("-subset-" + subset);
		
		return sw.toString();
	}
	
	public String toString() {
		List<String> lines = new ArrayList<String>();
		if (model != null)
			lines.add("model " + model.name);
		if (basket != null)
			lines.add("basket " + basket.name);
		if (property != null)
			lines.add("property " + property.getName());
		
		if (applier != null)
			lines.add("prediction for " + applier.compoundsProvider.getCompoundsNum() + " compounds");
		
		if (similarity != null && similarity > 0)
			lines.add("pairs similarity " + similarity);
		
		return StringUtils.join(lines, ", ");
	}
	
	public MMPStatsRequest(MMPSubset subset, double pValue) {
		this.subset = subset;
		this.pValue = pValue;
	}
	
	public MMPStatsRequest bootstrapStats() {
		bootstrapStatistics = true;
		return this;
	}
	
	public MMPStatsRequest setMinPairs(int minPairs) {
		this.minPairs = minPairs;
		return this;
	}
	
	public MMPStatsRequest setProperty(Property p)
	{
		this.property = p;
		return this;
	}
	
	public MMPStatsRequest() {
		
	}
}
