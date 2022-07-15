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

import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "validation-configuration")
public abstract class ValidationConfiguration extends ModelAbstractConfiguration {

	private static final long serialVersionUID = 1L;
	public static final int CLASSIC_BAGGING = 0;
	public static final int STRATIFIED = 2;

	public String taskName;
	public ModelAbstractConfiguration taskConfiguration;
	public int ensembleSize = 5;

	public Boolean keepIndividualPredictions;
	public Integer validationType;

	public static enum MixtureValidation {MIXTURE,  RECORD, ANIONS, CATIONS, MAXCOMP, COMPONENT, ALL};
	public MixtureValidation mixtureValidation; // by mixture for mixtures; if null no mixture validation is used

	@Override
	public String getInformativeName(){
		if(mixtureValidation == null) return "";
		switch(mixtureValidation) {
		case ANIONS:
			return " (anions)";
		case CATIONS:
			return " (cations)";
		case MAXCOMP:
			return " (maxcomp)";
		case COMPONENT:
			return " (component)";
		case ALL:
			return " (all)";
		case RECORD:
			return " (record)";
		default:
			return "";
		}
	}

	@Override
	public void setVersion(String version, boolean force)
	{
		super.setVersion(version, force);
		if(taskConfiguration != null)taskConfiguration.setVersion(version, force);
	}

	public int getBaggingType(){
		return validationType ==null?CLASSIC_BAGGING:validationType;
	}

	public boolean isRecordValidated() {
		return mixtureValidation != null && mixtureValidation == MixtureValidation.RECORD;
	}

	@Override
	public boolean isSupportRegression() {
		return true;
	}

	public boolean areMultiLearningData() {
		return taskConfiguration == null ?false:taskConfiguration.areMultiLearningData();
	}

	@Override
	public boolean isFeatureNet() {
		return taskConfiguration == null ?false:taskConfiguration.isFeatureNet();
	}

	@Override
	public boolean isForcedSingleTasklearning() {
		return taskConfiguration == null ?false:taskConfiguration.isForcedSingleTasklearning();
	}

	@Override
	public List<Integer> getOptions(){
		return taskConfiguration != null && taskConfiguration.getOptions() != null? taskConfiguration.getOptions() : optionsNumber;
	}

	@Override
	public void setTheOptions(List<Integer> options){
		taskConfiguration.setTheOptions(options);
	}

	@Override
	public ModelAbstractConfiguration getTrainingConfiguration() {
		return taskConfiguration;
	}

	public boolean compatibleWithDesalt(boolean isDesaulted) {
		if(mixtureValidation == null) return true;
		switch(mixtureValidation) {
		case ANIONS:
		case CATIONS:
		case MIXTURE:
		case MAXCOMP:
			return isDesaulted == false;
		default:
			return true;		
		}
	}

	public static void main(String[] args) throws Exception{
		ValidationConfiguration v = new CrossValidationConfiguration();
		v.mixtureValidation = MixtureValidation.MAXCOMP;
		System.out.println(v.getInformativeName());
	}

}
