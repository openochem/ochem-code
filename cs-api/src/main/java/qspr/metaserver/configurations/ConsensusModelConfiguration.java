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
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElementWrapper;
import javax.xml.bind.annotation.XmlRootElement;

import qspr.interfaces.ProvidedConditions;
import qspr.metaserver.util.ShortCondition;
import qspr.workflow.utils.QSPRConstants;

// Under construction / Midnighter on Jun 27, 2011
@XmlRootElement
public class ConsensusModelConfiguration extends ModelAbstractConfiguration implements ProvidedConditions
{
	private static final long serialVersionUID = 1L;
	public static enum ConsensusType {
		AVERAGE, OPTIMAL, BEST_MODEL, WEIGHTED_AVERAGE, RMSE_WEIGHTED
	}

	@XmlElementWrapper(name="individual-models")
	@XmlElement(name="individual-model")
	public List<IndividualModel> individualModels = new ArrayList<IndividualModel>();

	@XmlElement(name = "unit")
	public List<String> units;

	@XmlElementWrapper(name="models-for-properties")
	@XmlElement(name="models")
	public ArrayList<String> models;

	public List<Float[]> weights; // weights for RMSE-base average; size of weights == number of models

	@XmlElement(name = "classes")
	public List<Integer> classes; // weights for RMSE-base average; size of weights == number of models

	public ConsensusType type;

	public Boolean allowErrors;

	@XmlElementWrapper(name="external-descriptors")
	@XmlElement(name="external-descriptor")
	public List<ShortCondition> conditions;

	/**
	 * Child models should be added in the order of properties
	 * @param childModels
	 */
	public void addModelsForProperty(List<Long> childModels) {
		if(models == null) models = new ArrayList<String>(1);
		Collections.sort(childModels);
		String s = "";
		for(Long model:childModels)
			s += (s.length()>1?",":"") + model;
		models.add(s);
	}

	public void compact() {
		if(classes!= null) {
			boolean keep = false;
			for(Integer c:classes)
				if(c>2) { // at least one multiclass
					keep = true;
					break;
				}

			if(!keep)classes= null;			
		}
		if(models == null) return; // only if models per property
		List<HashSet<Long>> sets = getModelsForProperties();
		HashSet<Long> set = sets.get(0);
		for(HashSet<Long> s:sets)
			set.addAll(s);
		List<IndividualModel>  reducedSet = new ArrayList<IndividualModel>();		
		for(IndividualModel m:individualModels) 
			if(set.contains(m.id) && !reducedSet.contains(m))
				reducedSet.add(m);
		Collections.sort(reducedSet);
		individualModels = reducedSet;
		
		boolean ok = true;
		for(HashSet<Long> s:sets)
			for(IndividualModel m:individualModels)if(!s.contains(m.id)) {
				ok=false;
				break;
			}

		if(ok) { // if only one property or all models are identical
			models = null;
			type = ConsensusType.AVERAGE;
		}
		
	}

	public List<HashSet<Long>> getModelsForProperties() {
		if(models == null) return null;
		List<HashSet<Long>> set  = new ArrayList<HashSet<Long>>();
		for(String model:models) {
			HashSet<Long> s = new HashSet<Long>();
			for(String id:model.split(",")) {
				s.add(Long.parseLong(id));
			}
			set.add(s);
		}
		return set;
	}

	public void addModel(Long id)
	{
		for(IndividualModel model: individualModels)
			if((long)model.id == (long)id) return; // to avoid duplicates
		individualModels.add(new IndividualModel(id));
		Collections.sort(individualModels);
	}

	@Override
	public List<ShortCondition> getConditions() {
		return conditions;
	}

	@Override
	public boolean hasConditions() {
		return conditions != null;
	}

	@Override
	public boolean isTrainingConfiguration()
	{
		return false;
	}

	public String toString()
	{
		StringWriter buf = new StringWriter();
		buf.write(individualModels.size() + " individual models:\n");
		Collections.sort(individualModels); // to keep also old models sorted
		for (IndividualModel iModel : individualModels)
			buf.write("ModelID: " + iModel.id + "\n");

		buf.write(type+"\n");

		return buf.toString();
	}

	public static class IndividualModel implements Serializable, Comparable<Object>
	{
		private static final long serialVersionUID = 1L;

		@XmlAttribute
		public Long id;

		public IndividualModel(){
		}

		public IndividualModel(Long id){
			this.id=id;
		}

		@Override 
		public int compareTo(Object o) {
			IndividualModel f = (IndividualModel) o; 
			return (int)(this.id - f.id);
		}
	}

	@Override
	public boolean isSupportRegression() {
		return true;
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.CONSENSUS;
	}

	public boolean isModelSaved() {
		return models != null;
	}

}


