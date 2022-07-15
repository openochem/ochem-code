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

package qspr.metaserver.configurations;

import com.eadmet.utils.NumericalValueStandardizer;

abstract public class MultiLearningAbstractConfiguration extends ModelAbstractConfiguration implements Iterations{

	private static final long serialVersionUID = 1L;

	public Boolean noMultiLearning;

	public String experimentalParam;

	/**
	 * Global implicit values
	 */

	public Double[] implicitValues;

	/**
	 * The weighting schema of the training labels and classes. If empty, the equal weighting is assumed
	 */
	public LabelWeighting labelWeighting;

	@Override
	public boolean isForcedSingleTasklearning() {
		return (noMultiLearning != null && noMultiLearning) || isFeatureNet();
	}

	@Override
	public String toString()
	{
		return super.toString() + (isForcedSingleTasklearning() ?" STL":"") + (experimentalParam == null?"":" "+experimentalParam) + 
				(labelWeighting == null ? "": " " + labelWeighting.toString()) + (implicitValues == null? "": implicitValuesToString());
	}

	private String implicitValuesToString() {
		if(implicitValues == null) return "";
		String miss = "\nimplicitValues: ";
		for(Double val : implicitValues)
			miss += val == null? " -":" " +NumericalValueStandardizer.getSignificantDigits(val);
		return miss;
	}

	public Double [] getImplicit(){
		return implicitValues;
	}

	@Override
	public String getInformativeName() {
		return super.getInformativeName() + (getImplicit() != null ? " implicit":"") + " ";	
	}

}
