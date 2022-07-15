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

package qspr.modelling.configurations;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import qspr.interfaces.Descriptable;
import qspr.interfaces.ProvidedConditions;
import qspr.metaserver.configurations.DescriptorType;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsEmptyConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.NoDescriptors;
import qspr.metaserver.configurations.SelectionConfiguration;
import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.metaserver.util.ShortCondition;

@XmlRootElement(name = "cds-configuration")
public class CDSConfiguration implements Descriptable, ProvidedConditions
{
	private static final long serialVersionUID = 1L;
	public DescriptorsConfiguration descriptors = new DescriptorsConfiguration();
	public SelectionConfiguration selection = new SelectionConfiguration();
	//
	public StructureOptimisationConfiguration optimisationConfiguration;
	//
	@XmlElement(name = "externalDescriptors")
	public List<ExternalCondition> conditions;

	public ModelAbstractConfiguration modelConfiguration;

	// Some fields moved from here to Model object... to avoid serializing to database and stuff... NoS 2.12.09

	@Override
	public boolean hasConditions(){
		return conditions != null && conditions.size() > 0;
	}

	@SuppressWarnings({ "unchecked", "rawtypes" })
	@Override
	public List<ShortCondition> getConditions() {
		return hasConditions() ? (List<ShortCondition>) (List)conditions : null;
	}

	public boolean isCompatibleDescriptorsAndMethod(ModelAbstractConfiguration modelcfg) {
		if(hasDescriptors()) {
			if(modelcfg instanceof NoDescriptors) {
				if(!modelcfg.isSupportConditions()) return false;
				if(!onlyDescriptorsEmptyConfiguration()) return false;
			}
		}else
			if(!(modelcfg instanceof NoDescriptors)) return false;

		return true;
	}

	private boolean onlyDescriptorsEmptyConfiguration() {
		if(descriptors.types.size() == 0) return false;

		for(DescriptorType desc:descriptors.types)
			if(!(desc.configuration instanceof DescriptorsEmptyConfiguration))return false;
		return true;

	}

	public boolean hasDescriptors() {
		return descriptors.types.size() != 0;
	}

	public boolean hasMixtures() {
		return descriptors.mixtures != null;
	}

	public String toString()
	{
		String res = "";
		res += descriptors.toString()+"\n";
		if (conditions != null)
			res += conditions + "\n";
		res += selection+"\n";

		if (modelConfiguration != null)
			res += " " + modelConfiguration;

		if(optimisationConfiguration != null)
			res += " " + optimisationConfiguration;

		return res;
	}

	public Map<String, Object> getParameters() 
	{
		Map<String, Object> parameters = new HashMap<String, Object>();
		parameters.putAll(descriptors.getParameters());

		if (modelConfiguration instanceof Descriptable)
			parameters.putAll(((Descriptable) modelConfiguration).getParameters());

		return parameters;
	}

}
