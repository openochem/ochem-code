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

import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.FileUtils;

import qspr.metaserver.util.DataScaling;
import qspr.util.ClassCompressor;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

/**
 * Represents the configuration of a classical predictive model that: - takes input samples (descriptors) and labels (experimental values) - gives back
 * predictions - gives back a predictive model
 */
abstract public class ModelAbstractConfiguration implements Serializable {
	private static final long serialVersionUID = 1L;

	public Boolean featureNet; 

	public Boolean saveModels;

	public DataTable iterations;

	public enum DataType {
		REGRESSION, CLASSIFICATION, MULTIMIXED, MULTICLASSIFICATION, MULTIREGRESSION
	};

	/**
	 * Store your saved model as this object
	 */
	public CompressedObject<Object> savedmodel; // since 10.08 this is new way to store models
	public CompressedObject<DataScaling> scaleX, scaleY; // since 26.08 this is the way to store globalnormalisation

	public CompressedObject<Object> uploadedModelData;

	public ScalingType scaleTypeX, scaleTypeY; // scaling types for x and y values

	public List<Integer> optionsNumber; // Size = number of modeled properties, optionsNumber[n] = number of classes for property "n", 1 for regression; enumeration is done during model preparation
	public Map<Integer, Integer> qualitativeDescriptors; // Column Number -> Max ValueID present in the training set
	public HashMap<String, CompressedObject<Serializable>> dmModels;

	public Integer seed;

	public PredictionScenario predictionScenario;

	public CompressedObject<Object> featureNetConfig;

	public ModelAbstractConfiguration() {
		super();
	}

	public int getSeed() {
		return seed == null?QSPRConstants.SEED:seed;
	}

	public String versionOCHEM = null;

	private Boolean skipSizeCheck;

	public void setVersion(String version, boolean force){
		if(versionOCHEM == null || force)
			versionOCHEM = version != null? version.replaceAll("\\r\\n", "").trim():null;
	}

	public List<Integer> getOptions(){
		return optionsNumber;
	}

	/**
	 * To avoid using setOptions, which will be considered as get/set
	 * @param options
	 */
	public void setTheOptions(List<Integer> options){
		optionsNumber = options;
	}

	public String toString() {
		String s = versionOCHEM == null || versionOCHEM.startsWith("v.3.")?"":versionOCHEM + " ";
		s += scaleTypeX != null && scaleTypeX != ScalingType.NONE ? "\nscale X: " + scaleTypeX : "";
		s += scaleTypeY != null && scaleTypeY != ScalingType.NONE ? "\nscale Y: " + scaleTypeY : "";
		if(isFeatureNet()) s += " FeatureNet";
		return s;
	}

	public void cleanBulkyStuff() {
		savedmodel = null;
		scaleX = scaleY = null;
		dmModels = null;
		featureNetConfig = null;
	}

	public Object getSavedModelAsObject() {
		return savedmodel.get();
	}

	public byte[] getSavedModelAsBytes() {
		return (byte[])getSavedModelAsObject();
	}

	/**
	 * 
	 * @param scale
	 * @param type
	 *            0 is for x, 1 is for y
	 */

	public void addScaleY(DataScaling scale) {
		scaleY = new CompressedObject<DataScaling>();
		scaleY.set(scale);
	}

	public void addScaleX(DataScaling scale) {
		scaleX = new CompressedObject<DataScaling>();
		scaleX.set(scale);
	}

	public boolean saveModels() {
		return saveModels == null || saveModels;
	}

	public boolean isModelSaved() {
		if(isFeatureNet())return featureNetConfig != null && savedmodel != null;
		return savedmodel != null;
	}

	public boolean isTrainingConfiguration() {
		return savedmodel == null;
	}

	public boolean isUploadedConfiguration() {
		return (uploadedModelData != null);
	}

	final public void storeModel(File f ) throws IOException {
		if(!saveModels())return;
		if(f.exists() && f.isDirectory()) {
			String s = FileUtils.createTarFile(f.getCanonicalPath());
			if(s.length() > QSPRConstants.MAXMODELSIZE) saveModels = false;
			else
			{
				setModel(FileUtils.getFileAsBytes(s));
				org.apache.commons.io.FileUtils.deleteDirectory(f);
			}
			return;
		}

		if(!f.exists() || f.length() == 0) throw new IOException("Model file is absent: " + f.getCanonicalPath());
		if(f.length() > QSPRConstants.MAXMODELSIZE)saveModels = false; // at least we can try to apply model 
		else{
			setModel(FileUtils.getFileAsBytes(f.getCanonicalPath()));
			f.delete();
		}
	}

	public void setModel(Object themodel) {
		if (themodel == null)
			savedmodel = null;
		else {
			savedmodel = new CompressedObject<Object>();
			savedmodel.set(themodel);
		}
	}


	/**
	 * Return list of descriptors that were selected Methods which do variable selection should override this method
	 * 
	 * @return
	 */

	public List<Integer> getSelectedDescriptors() {
		return null;
	}

	public ModelAbstractConfiguration getDeepCopy() {
		// an unoptimal, but universal way of copying a class - serialize and unserialize / (c) Midnighter
		return (ModelAbstractConfiguration) ClassCompressor.byteToObject(ClassCompressor.objectToByte(this));
	}

	public enum ScalingType {
		NONE, STANDARDIZE, RANGE, RANGE_MINUS1_PLUS1;
	}

	public PredictionScenario getPredictionScenario() {
		return predictionScenario == null ? PredictionScenario.PREDICTION_ONLY : predictionScenario;
	}

	public void addDM(String dm) {
		addDM(dm, null);
	}

	public void addDM(String dm, Serializable dmModel) {
		if (dmModels == null)
			dmModels = new HashMap<String, CompressedObject<Serializable>>();

		if (dmModel == null)
			dmModels.put(dm + "#" + dm.length(), null);
		else
			dmModels.put(dm + "#" + dm.length(), new CompressedObject<Serializable>(dmModel));
	}

	public ModelAbstractConfiguration getOptimisedConfiguration() throws IOException {
		return this;
	}

	public enum PredictionScenario {
		/**
		 *  Only predictions, even for training set compounds
		 */
		PREDICTION_ONLY(0), 

		/**
		 * CV predictions for training set + real predictions for new compounds
		 */
		CROSS_VALIDATED_PREDICTIONS(1),

		/**
		 *  Experimental values for training set compounds
		 */
		PREDICTION_AND_TABLE(2),

		/**
		 *  use to update ASNN/Bagging predictions to find nearest neighbors; internal option which is used only inside of these servers
		 */
		DISTANCE_ONLY(3);

		private byte flag;

		private PredictionScenario(int flag)
		{
			this.flag = (byte) flag;
		}

		public byte getValue()
		{
			return (byte) flag;
		}

		static PredictionScenario fromValue(byte arg)
		{
			switch (arg)
			{
			case 0:
				return PREDICTION_ONLY;
			case 1:
				return CROSS_VALIDATED_PREDICTIONS;
			case 2:
				return PREDICTION_AND_TABLE;
			case 3:
				return DISTANCE_ONLY;
			default:
				throw new RuntimeException("Invalid argument for PredictionScenario " + arg);
			}
		}

	}

	public boolean isSupportRegression(){
		return false;
	}

	public boolean isForcedSingleTasklearning() {
		return false;
	}

	public boolean isSupportConditions(){
		return isSupportDescriptors();
	}

	public boolean isSupportDescriptors() {
		return true;
	}

	public ModelAbstractConfiguration.DataType getExtendedDataType() {
		Boolean regression = false, classification = false, multiclass = false;

		if (optionsNumber == null) return DataType.REGRESSION; // required for testing when options are not enumerated

		for (int i = 0; i < optionsNumber.size(); i++)
			if (optionsNumber.get(i) > 1) {
				classification = true;
				if (optionsNumber.get(i) > 2)
					multiclass = true;  // such cases are considered as multilearning
			}
			else
				regression = true;

		if(classification  && regression) return DataType.MULTIMIXED; // multiple class classifications are treated as MIXED types requiring multilearning
		if(regression && optionsNumber.size() > 1) return DataType.MULTIREGRESSION;
		if(multiclass || optionsNumber.size() > 1) return DataType.MULTICLASSIFICATION;
		if(classification) return DataType.CLASSIFICATION;
		return DataType.REGRESSION;
	}

	public Integer OutputValues() {

		if (optionsNumber == null)throw new UserFriendlyException("Something is wrong: optionsNumber are not enumerated!");

		int count = 0;
		for (int i = 0; i < optionsNumber.size(); i++)
			count += optionsNumber.get(i);

		return count;
	}

	public boolean containAlsoRegressionData() {
		if (optionsNumber == null) return true;
		for (int i = 0; i < optionsNumber.size(); i++)
			if (optionsNumber.get(i) == 1) return true;
		return false;
	}

	public boolean containMultiClassData() {
		if (optionsNumber == null) return false;
		for (int i = 0; i < optionsNumber.size(); i++)
			if (optionsNumber.get(i) > 2)
				return true;
		return false;
	}

	public boolean areMultiLearningData() {
		if (optionsNumber == null) return false;
		return optionsNumber.size() > 1;
	}

	public boolean areClassificationData() {
		if (optionsNumber == null)throw new UserFriendlyException("Something is wrong. optionsNumber are not enumerated!");
		for (int i = 0; i < optionsNumber.size(); i++)
			if (optionsNumber.get(i) > 1)
				return true;
		return false;
	}

	abstract public String getDefaultName();

	private String getSTL() {
		return (isFeatureNet()?" (FN)":isForcedSingleTasklearning()?" (STL)":"");
	}

	public String getInformativeName() {
		return getDefaultName() + getSTL() + (saveModels()?"":QSPRConstants.UNSAVED)+ " ";
	}

	public boolean isSupportPredicates() {
		return false;
	}

	public ModelAbstractConfiguration getTrainingConfiguration() {
		return this;
	}

	public int requireMinimumRecords() {
		return 10;
	}

	public boolean isLarge(){
		return false;
	}

	public boolean isFeatureNet() {
		return featureNet != null && featureNet && areMultiLearningData();
	}

	public boolean skipSizeCheck() {
		return isForcedSingleTasklearning() || (skipSizeCheck != null && skipSizeCheck);
	}

	public void setFeatureNet(Boolean b){
		featureNet = b;
	}
}
