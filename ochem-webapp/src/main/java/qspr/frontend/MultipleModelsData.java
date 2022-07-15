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

package qspr.frontend;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import com.eadmet.utils.NumericalValueStandardizer;

import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.Property;
import qspr.modelling.SetStatistics;
import qspr.workflow.utils.QSPRConstants;

/**
 * The "view bean" storing the statistics about multiple models
 * @author midnighter
 *
 */
@XmlRootElement(name = "multiple-models")
public class MultipleModelsData
{
	public Basket trainingSet;

	public List<String> methodCodes = new ArrayList<String>();
	public List<String> descriptorCodes = new ArrayList<String>();
	public transient List<String> methodHashes = new ArrayList<String>();
	public transient List<String> descriptorHashes = new ArrayList<String>();
	public transient List<String> trainingSetHashes = new ArrayList<String>();

	public String validation;
	public String validationTitle;
	public boolean hasAnyModels = false;

	final static int COUNTS = 8; // filed for COUNTS


	public ModelData[][] modelData ;

	String order(String desc) {
		return desc.startsWith("SMILES")?desc.toLowerCase():desc;
	}

	public void sortDescriptors() {
		if(modelData == null)return;

		TreeSet<String> sortedDesc = new TreeSet<String>();
		for(String desc:descriptorCodes)
			sortedDesc.add(order(desc));
		Map<Integer,Integer> newToOld= new HashMap<Integer,Integer>();
		int n=0;
		for(String desc:sortedDesc){
			TreeSet<String> sortedHashes =  new TreeSet<String>();
			for(int j=0;j<descriptorCodes.size();j++)
				if(order(descriptorCodes.get(j)).equals(desc))sortedHashes.add(descriptorHashes.get(j));
			for(String hash:sortedHashes)
				for(int j=0;j<descriptorHashes.size();j++)
					if(descriptorHashes.get(j).equals(hash))newToOld.put(n++,j);
		}

		List<String> descriptorCodesNew = new ArrayList<String>();
		List<String> descriptorHashesNew = new ArrayList<String>();
		ModelData[][] modelDataNew = new ModelData[modelData.length][modelData[0].length];

		for(int i=0;i<modelData.length;i++)
			for (Map.Entry<Integer, Integer> pair : newToOld.entrySet()) {
				modelDataNew[i][pair.getKey()] = modelData[i][pair.getValue()];
			}

		for (Map.Entry<Integer, Integer> pair : newToOld.entrySet()) {
			descriptorCodesNew.add(descriptorCodes.get(pair.getValue()));
			descriptorHashesNew.add(descriptorHashes.get(pair.getValue()));
		}

		modelData = modelDataNew;
		descriptorCodes = descriptorCodesNew;
		descriptorHashes = descriptorHashesNew;
	}

	public void addMethod() {
		int size1 = modelData == null ? 1 : modelData.length+1;
		int size2 = modelData == null ? 1 : modelData[0].length;
		copyAll(size1, size2);
	}

	public void addDescriptors() {
		int size1 = modelData == null ? 1 : modelData.length;
		int size2 = modelData == null ? 1 : modelData[0].length + 1;
		copyAll(size1, size2);
	}

	void copyAll(int size1, int size2) {

		ModelData[][] modelDat = new ModelData[size1][size2];

		if(modelData != null)
			for(int i = 0; i<modelData.length;i++)
				for(int j = 0; j< modelData[0].length;j++)modelDat[i][j] = modelData[i][j];
		modelData = modelDat;
	}

	@XmlElement(name = "models-count")
	public int getModelCount()
	{
		if(modelData == null) return 0;
		int count = 0;
		for (int i = 0; i < modelData.length; i++)
			for (int k = 0; k < modelData[i].length; k++)
				if (modelData[i][k] != null)
					count += modelData[i][k].count;
		return count;
	}

	@XmlElement
	protected int getTrainingSetVersionCount() {
		return trainingSetHashes.size();
	}

	/**
	 * Return type of model(s) 
	 *  0 - regression
	 *  1 - classification
	 *  2 - mixed
	 * @return
	 */

	@XmlElement
	public ModelType getModelType(){ // determines type based on the first property does not correctly work for exporting of multiple model statistics

		Boolean classification = null;

		for(Property prop: trainingSet.getProperty()) {
			if(classification != null && classification != prop.isQualitative()) return ModelType.mixed;
			classification = prop.isQualitative();
		}

		return classification ? ModelType.classification : ModelType.regression;
	}

	public static enum ModelType{
		regression, classification, mixed
	}


	public static class ModelData
	{
		public STATUS status;
		public String detailedStatus;

		public static enum STATUS {
			published, saved, init, running, ready, deleted, invalid, error, // used statuses
			kill, assigned, in_queue; // replaced statuses
		};


		public Long modelId;
		public Long publicModelId;
		public Long pendingTaskId;
		public String modelName;
		public String description;
		public int count = 1;
		public List<ModelData> moreModels;
		public List<SetStats> stats = new ArrayList<SetStats>();


		public void addAverage() {
			Map<String,Integer> baskets = new HashMap<String, Integer>();

			for(SetStats set: stats)
				if(baskets.containsKey(set.id))baskets.put(set.id, baskets.get(set.id) + 1);
				else
					baskets.put(set.id, 1);

			for(String basket: baskets.keySet()){  // for each basket

				if(baskets.get(basket) == 1) continue;

				Double vals[][] = new Double[baskets.get(basket)][];

				int i = 0;
				SetStats found = null;
				for(SetStats set: stats)if(basket.equals(set.id)) try{ 
					vals[i++] = set.toValues();
					found = set;
				}catch(Exception e) {

				}

				int length  = 0;
				for(i = 0; i < vals.length; i++)
					if(vals[i] != null) {
						length = vals[i].length;
						break;
					}

				System.out.println(basket +" " + vals.length+ " len: " + length);
				
				if(length == 0) return;

				Double average[] = new Double [length], count [] =  new Double [length];

				for(i = 0; i < vals.length; i++)
					if(vals[i] != null)
						for(int j = 0; j < length; j++)if(vals[i][j] != null && Double.isFinite(vals[i][j])){
							if(average[j] == null || vals[i][COUNTS] == null || vals[i][COUNTS] == 0.) average[j] = count[j] = 0.;
							average[j] += vals[i][j];
							count[j]++;
						}

				for(int j = 0; j < COUNTS; j++)if(count[j] != null) // after COUNTS no averaging is required
					average[j] /= count[j];
				else
					average[j] = null;

				SetStats aver = new SetStats(average);
				aver.basketId = found.basketId;
				aver.id = found.id;

				stats.add(aver);
			}

		}

		public void setStatus(String theStatus){

			try{
				status = STATUS.valueOf(theStatus.replace(' ', '_'));
			}catch(Exception e){
				System.out.print("Exception: unknown status "+theStatus);
				status = STATUS.invalid;
			}

			switch(status){
			case kill: status = STATUS.error; break;
			case assigned: status = STATUS.running; break;
			case in_queue: status = STATUS.init; break;
			default: break;
			}
		}

		public void setModel(Model model)
		{
			modelId = model.id;
			modelName = model.name + " " + model.publicId;
			description = model.description;
			publicModelId = model.publicId;
		}

		public void addStatistics(SetStatistics ss, Model m, boolean first)
		{
			try
			{
				stats.add(new SetStats(ss, m, first));
			} catch (Exception e)
			{
				e.printStackTrace();
				status = ModelData.STATUS.invalid;
				detailedStatus = e.getMessage();
			}
		}

		public SetStats getStatsByBasketID(String basketID) {
			SetStats found = null; // export average or last value
			for (SetStats ss : stats)
				if (ss.getBasketSetNameId().equals(basketID))
				{
					found = ss; 
					break;
				}
			return found;
		}

		public void addAnotherModel(ModelData mData)
		{
			if (moreModels == null)
				moreModels = new ArrayList<MultipleModelsData.ModelData>();
			moreModels.add(mData);
			count++;
		}
	}

	public static class SetStats
	{
		public String rmse;
		public String r2;
		public String q2;
		public String mae;
		public String average;
		public String accuracy;
		public String balancedAccuracy;
		public String auc;
		public String size;
		public String records;
		public String errors;

		@XmlAttribute
		public String id;

		public Long basketId;

		public SetStats()
		{

		}

		public String getBasketSetNameId() {
			return "" + basketId + "_" + id;
		}

		/**
		 * Uses values which are provided for each property just the function below: Double[] toValues()
		 * @param vals
		 */
		public SetStats(Double [] vals)
		{
			if(vals[0] != null) rmse = "(" + NumericalValueStandardizer.getSignificantDigits(vals[0]) + ")";
			if(vals[1] != null) r2 = "(" + NumericalValueStandardizer.getSignificantDigits(vals[1]) + ")";
			if(vals[2] != null) q2 = "(" + NumericalValueStandardizer.getSignificantDigits(vals[2]) + ")";
			if(vals[3] != null) mae = "(" + NumericalValueStandardizer.getSignificantDigits(vals[3]) + ")";
			if(vals[4] != null) average = "(" + NumericalValueStandardizer.getSignificantDigits(vals[4]) + ")";
			if(vals[5] != null) accuracy = "(" + NumericalValueStandardizer.getSignificantDigits(vals[5]) + ")";
			if(vals[6] != null) balancedAccuracy = "(" + NumericalValueStandardizer.getSignificantDigits(vals[6]) + ")";
			if(vals[7] != null) auc = "(" + NumericalValueStandardizer.getSignificantDigits(vals[7]) + ")";
			if(vals[COUNTS] != null) records = "(" + Math.round(vals[COUNTS]) + ")";
			if(vals[9] != null) errors = "(" + Math.round(vals[9]) + ")";
		}

		Double[] toValues() {
			Double vals[] = new Double [10];
			if(rmse != null )vals[0] = Double.valueOf(rmse);
			if(r2 != null )vals[1] = Double.valueOf(r2);
			if(q2 != null )vals[2] = Double.valueOf(q2);
			if(mae != null )vals[3] = Double.valueOf(mae);
			if(average != null )vals[4] = Double.valueOf(average);
			if(accuracy != null )vals[5] = Double.valueOf(accuracy);
			if(balancedAccuracy != null )vals[6] = Double.valueOf(balancedAccuracy);
			if(auc != null )vals[7] = Double.valueOf(auc);
			if(records != null )vals[COUNTS] = Double.valueOf(records);
			if(errors != null )vals[9] = Double.valueOf(errors);
			return vals;
		}


		public SetStats(SetStatistics ss, Model m, boolean first)
		{
			int n = ss.points.size();

			rmse = n == 0 ?"-":ss.getRmseStr();
			mae = n == 0 ?"-":ss.getMaeStr();
			average = n == 0 ?"-":ss.getAverageStr();
			r2 = n == 0 ?"-":ss.getR2Str();
			q2 = n == 0 ?"-":ss.getQ2Str();
			records = "" + ss.points.size();
			errors = "" +ss.errorCount;
			if (ss.getAccuracyStr() != null)
			{
				accuracy = n == 0 ?"-":ss.getAccuracyStr();
				balancedAccuracy = n == 0 ?"-":ss.getBalancedAccuracyStr();
				auc = n == 0 ?"-":ss.getAucAccuracyStr();
				r2 = n == 0 ?"-":ss.getMCCStr();
				rmse = mae = average = q2 = null;
			}

			id = ss.setId;
			basketId = ss.basketId;

			int mB = 1024 * 1024;
			size = first ? (m.size > mB ? NumericalValueStandardizer.getSignificantDigitsStr(m.size / mB, 0) + "Mb" : NumericalValueStandardizer.getSignificantDigitsStr(m.size / 1024, 0)  + "Kb") : "";
		}
	}

	public MultipleModelsData setValidation(String validation, String title)
	{
		this.validation = validation;
		this.validationTitle = title;
		return this;
	}

	private void cleanWhole(List<SetStats> a){
		for(int j=0;j<a.size();j++)
			if(a.get(j).id.equals(QSPRConstants.WHOLE))
				a.remove(j);
	}

	public void cleanWhole() {
		if(modelData != null)for (int i = 0; i < modelData.length; i++)
			if(modelData[i] != null)for (int k = 0; k < modelData[i].length; k++) {
				if (modelData[i][k] != null){
					cleanWhole(modelData[i][k].stats);
					if (modelData[i][k].moreModels != null)
						for(ModelData dat: modelData[i][k].moreModels)
							cleanWhole(dat.stats);
				}
			}
	}

	public void addAverage() {
		if(modelData != null)
			for(ModelData [] mode : modelData) 
				if(mode != null) for(ModelData m : mode) 
					if(m != null) {
						m.addAverage();
						if(m.moreModels != null)
							for(ModelData mm: m.moreModels)mm.addAverage();
					}
	}

	/**
	 * All validation sets should get the same label
	 */

	private void anonymyseValidationSets(List<SetStats> stats) {
		for (SetStats ss : stats) {
			if (ss.id.startsWith(QSPRConstants.VALIDATION))
				ss.id = QSPRConstants.VALIDATION; // "thus, validation0" -> QSPRConstants.VALIDATION. All validation sets are same for us here
		}
	}

	public void anonymyseValidationSets() {
		for(int i = 0; modelData != null && i<modelData.length;i++)
			for(int j = 0; modelData[0] != null && j< modelData[0].length;j++)
			{
				ModelData mData = modelData[i][j];
				if(mData != null) {
					anonymyseValidationSets(mData.stats);
					if (mData.moreModels != null)
						for(ModelData dat: mData.moreModels)
							anonymyseValidationSets(dat.stats);
				}
			}
	}

}
