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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

import qspr.util.ClassCompressor;
import qspr.workflow.datatypes.StringList;

@XmlRootElement(name = "selection-configuration")
public class SelectionConfiguration implements Serializable 
{
	private static final long serialVersionUID = 1L;

	public StringList descriptors;
	public StringList fixedDescriptors;

	/*
	 * Default settings are as in browser
	 */

	public int numDifferentValues = 2;
	public double correlationThreshold = 0.95;
	public boolean useUFS = false;
	public Boolean useAUTO;
	public CompressedObject<Object> autoEncoder;
	public double stdThreshold = 0.01;
	public int maximumValueThreshold = 999999;

	public String toStringIgnoreSelected() {
		String res = "Correl. limit: " + correlationThreshold;
		if (numDifferentValues != 2)
			res += " Unique values: " + numDifferentValues + ", ";
		res += " Variance threshold: "+stdThreshold + ", ";
		res += " Maximum value: "+maximumValueThreshold + ", ";
		if (useUFS)
			res += "using UFS";

		return res;
	}

	public String toString() {
		String res = "Correl. limit: " + correlationThreshold;
		if (numDifferentValues != 2)
			res += " Unique values: " + numDifferentValues + ", ";
		res += " Variance threshold: "+stdThreshold + ", ";
		res += " Maximum value: "+maximumValueThreshold + ", ";
		if (useUFS)
			res += "using UFS";

		if (descriptors != null && descriptors.values.size() > 0)
			res += ", " + descriptors.values;

		return res;
	}

	/*
	 * Writes cfg file and also checks whether with these parameters we should
	 * indeed use filtering
	 */

	public boolean doFiltering() 
	{
		if (stdThreshold > 0)
			return true;
		if (maximumValueThreshold < 999999)
			return true;
		if (numDifferentValues > 1)
			return true;
		if (correlationThreshold > 0 && correlationThreshold < 1)
			return true;
		return false;
	}

	public void setNoFiltering() {
		stdThreshold = 0;
		maximumValueThreshold = Integer.MAX_VALUE;
		numDifferentValues = 0;
		correlationThreshold = 0;
	}

	public int getDescriptorsSize(){
		return descriptors == null || descriptors.values == null? 0 : descriptors.values.size();
	}

	public List<String> descriptorAsStrings() {
		return getDescriptorsSize() == 0? new ArrayList<String>() : descriptors.values;
	}

	/**
	 * Add descriptors and trims them before adding to the list
	 * @param theDescriptors
	 */

	public void storeDescriptorAsStrings(List<String> theDescriptors) {
		descriptors = new StringList();
		for(String descriptor: theDescriptors){
			descriptors.values.add(descriptor.trim());
		}

	}

	public SelectionConfiguration getDeepCopy() {
		return (SelectionConfiguration) ClassCompressor.byteToObject(ClassCompressor.objectToByte(this));
	}

}
