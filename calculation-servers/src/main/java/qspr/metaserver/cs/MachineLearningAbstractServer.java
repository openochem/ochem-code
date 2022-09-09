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

package qspr.metaserver.cs;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import qspr.interfaces.DataDrivenConfiguration;
import qspr.metaserver.configurations.BaggingConfiguration;
import qspr.metaserver.configurations.CompressedObject;
import qspr.metaserver.configurations.Iterations;
import qspr.metaserver.configurations.LSSVMGConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.configurations.NoDescriptors;
import qspr.metaserver.configurations.ValidationConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.protocol.Task.DebugLevel;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.CalculationTask;
import qspr.workflow.utils.CalculationTaskSet;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.utils.NumericalValueStandardizer;
import com.eadmet.utils.OCHEMUtils;

public abstract class MachineLearningAbstractServer extends WorkflowNodeServer
{
	private static final String[] DEFAULT_MESSAGES = { "EPOCH:","TRAIN SCORE:","VALIDATION SCORE:"};
	protected int batchApplySize = Integer.MAX_VALUE/2; // maxmimum apply size of a model
	final int MIN_BATCH_SIZE= 10000;

	protected DataTable iterations;

	private String messages[] = null;

	private double shift = 0.;

	// The original data are stored to have a possibility to identify sizes
	public boolean allowLocalCalculations = false;

	/**
	 * Train a predictive model
	 * 
	 * @param dtDescriptors - descriptors for both the training and test sets
	 * @param dtExpValues - labels for the training set
	 * @param configuration - configuration for training the model and for storing the newly created model
	 * @return DataTable of predictions for the training and test sets
	 * @throws Exception
	 */
	abstract protected DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws Exception;

	/**
	 * Apply a predictive model to a subset of molecules
	 * 
	 * @param dtDescriptors descriptors for the predicted compounds
	 * @param receivedConf
	 * @return
	 * @throws Exception
	 */
	abstract protected DataTable applyModelBatch(DescriptorsTable dtDescriptors, ModelAbstractConfiguration receivedConf) throws Exception;

	void initIterations() {

		iterations = new DataTable();
		iterations.addRow();
		for(int i=0;i<messages.length;i++) {
			messages[i]=messages[i].toUpperCase();
			iterations.addColumn(DEFAULT_MESSAGES[i]);
			iterations.setValue(0,i,0.); // always to set to have some values
		}
	}

	@Override
	public void scanStatus(String newMessage) 
	{
		if(messages == null) {
			messages = getMessages();
			if(messages.length != 3 )throw new CriticalException("Number of messages to monitor should be 3");
		}

		newMessage = newMessage.toUpperCase();

		if(iterations == null)
			initIterations();

		boolean newEpoch = false;

		for(int i =0;i<messages.length;i++) 
			if(newMessage.contains(messages[i])){
				iterations.setValue(i, parseStatus(newMessage,messages[i]));
				if(iterations.getRowsSize() == 3) 
					for(int c=1;c<iterations.columns.size();c++) // filling in missed values
						if(iterations.getValue(0, c) == null)
							iterations.setValue(0, c,iterations.getValue(1, c));
				if(i == 0) newEpoch = true;
			}

		if(newEpoch) {
			int rows = iterations.getRowsSize();

			String m = "";
			for(int l =0; l<DEFAULT_MESSAGES.length && l <iterations.columns.size(); l++)
				if(iterations.getValue(l) != null)
					m += " " + DEFAULT_MESSAGES[l] + " " + iterations.getValue(l);
			super.setStatus(m);
			if(rows > 2) {
				double now = (double) iterations.getValue(0);
				double before = (double) iterations.getValue(rows-2, 0);

				if(before > now && shift == 0) // fix for restarts of the transformer
					shift = before;
				iterations.setValue(now + shift);
			}

			iterations.addRow();
			for(int c=0;c<iterations.columns.size();c++)
				iterations.setValue(c,iterations.getValue(rows-1,c)); // always to set to have some values
		}

	}

	private double parseStatus(String message, String keyword) {
		try {
			message = message.substring(message.indexOf(keyword)+keyword.length());
			message = message.replaceAll("[^0-9\\.\\+\\-]"," ").trim();
			String m[] = message.split("\\s+");
			//super.setStatus("**** " + keyword + " " + m[0]+ " in:  " + message);
			return Double.parseDouble(NumericalValueStandardizer.getSignificantDigits(Double.parseDouble(m[0])));
		}catch(Exception ee) {
			return 0.;
		}
	}

	/**
	 * @return messages that are processed to create  performance table
	 */
	protected String[] getMessages() {
		return DEFAULT_MESSAGES;
	}

	final DataTable applyModel(DescriptorsTable dtDescriptors, ModelAbstractConfiguration receivedConf) throws Exception {

		String name = "Applying " + (supportedTaskType.equalsIgnoreCase(receivedConf.getInformativeName())?supportedTaskType:supportedTaskType + "/" + receivedConf.getInformativeName()) + " to ";

		setStatus(name + dtDescriptors.getDataSize() + " molecules");

		int batchApplySize = this.batchApplySize;

		if(receivedConf instanceof NoDescriptors) {
			if(((NoDescriptors)receivedConf).getAugmentApply() >1)
				batchApplySize = batchApplySize/(((NoDescriptors)receivedConf).getAugmentApply());
			batchApplySize = batchApplySize < MIN_BATCH_SIZE ? MIN_BATCH_SIZE : batchApplySize;

			batchApplySize = batchApplySize > this.batchApplySize ? this.batchApplySize : batchApplySize;
		}

		int compressing [] = new int[dtDescriptors.getDataSize()]; 

		dtDescriptors = dtDescriptors.compressDescriptorsForApply(compressing);

		DataTable res =  fixModelBatch(dtDescriptors, receivedConf, 2*batchApplySize);
		if(compressing.length != res.getRowsSize()) // if there was a compressing
			res = res.decompress(compressing);

		return res;
	}

	DataTable fixModelBatch(DescriptorsTable dtDescriptors, ModelAbstractConfiguration receivedConf, int batchSize) throws Exception {

		DataTable res = null;
		if(batchSize >= dtDescriptors.getDataSize())
			try {
				return applyModelBatch(dtDescriptors, receivedConf);
			}catch(Exception e) {
				System.out.println(e.getMessage());
				if(dtDescriptors.getDataSize() == 1) {
					res = new DataTable(true);
					res.addRow().setError(e.getMessage());
					return res;
				}
				if(debug == DebugLevel.NONE) // not further propagation for non debug level
					throw e;
			}

		long start = Calendar.getInstance().getTimeInMillis();

		// more than 1 row is still available

		batchSize = dtDescriptors.getDataSize()/2 < batchSize/2? dtDescriptors.getDataSize()/2: batchSize/2 ;
		for(int i = 0; i < dtDescriptors.getDataSize() ; i += batchSize) {
			setStatus("Applying batch = " + (i/batchSize +1) + " out of n = " + (dtDescriptors.getDataSize()/batchSize + 1) +  " with batchSize = " + batchSize + " elapsed = "+ (Calendar.getInstance().getTimeInMillis()-start)/1000 + " sec");

			long size = i + batchSize < dtDescriptors.getDataSize() ? i + batchSize : dtDescriptors.getDataSize();
			DescriptorsTable mols = dtDescriptors.getSliceDescriptorsTable(i, (int)size);
			DataTable resnew = fixModelBatch(mols, receivedConf, batchSize);

			if(res == null) res = resnew;
			else
				res.addRowsFrom(resnew);

		}

		return res;

	}


	protected DataTable processUploadedModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration configuration)
			throws Exception
	{
		throw new CriticalException("Unfortunately, this machine learning method doesn't support model upload.");
	}


	void setIterations(ModelAbstractConfiguration receivedConf) {
		int iter = -1;

		for(int i =0; receivedConf.iterations != null && i < receivedConf.iterations.getRowsSize();i++) {
			Double val;
			if(( val = (Double)receivedConf.iterations.getValue(i, 0)) != null) {
				if(iter < val)iter = val.intValue();
			}
		}
		if(receivedConf instanceof Iterations && iter >= 0) { 
			setStatus("Used iterations = iter");
			((Iterations)receivedConf).setIterations(iter);
		}

	}

	@Override
	final public WorkflowNodeData calculate(WorkflowNodeData wndInput, Serializable receivedConf) throws Exception
	{
		ModelAbstractConfiguration configuration = (ModelAbstractConfiguration) receivedConf;
		init();

		DataTable dtDescriptors = wndInput.ports.get(0);
		DataTable dtValues = wndInput.ports.size() > 1 ? wndInput.ports.get(1) : null;

		if(dtDescriptors.getRowsSize() != dtDescriptors.getRowsNoErrorsSize())
			throw new CriticalException("The descriptor table contains descriptors with all errors. Try diffrent descriptors or/and verify your data.");	

		boolean training = configuration.isTrainingConfiguration();

		DataTable dtResult  = training? calculatePreFillTrain(dtDescriptors, dtValues, configuration): calculateApply(dtDescriptors, configuration);

		calculateDMsSafely(wndInput, dtResult, configuration, training);

		WorkflowNodeData wndResult = new WorkflowNodeData(dtResult);

		if (training)
		{
			configuration = configuration.getOptimisedConfiguration(); // since for CV we need to get taskConfiguration here!

			if (!configuration.saveModels())
				configuration.cleanBulkyStuff();

			wndResult.addPort(new DataTable(configuration));
		}

		if(dtResult.getRowsNoErrorsSize() != 0) // at least something was calculated
			cleanup();

		for(int i=0;i<dtDescriptors.getRowsSize();i++)
			dtDescriptors.getRow(i).attachments.remove(QSPRConstants.SDF_COLUMN);

		return wndResult;
	}


	final private DataTable calculateApply(DataTable dtDescriptors, ModelAbstractConfiguration configuration) throws Exception
	{
		if(dtDescriptors.getRowsSize() == 0)
			return new DataTable();

		/**
		 * All rows will contain information whether they are training or have validation type 
		 */
		for (int i = 0; i < dtDescriptors.getRowsSize(); i++)
		{
			AbstractDataRow r = dtDescriptors.getRow(i);
			r.removeAttachment(QSPRConstants.VALIDATION);
			r.addAttachment(QSPRConstants.VALIDATION, true);
		}

		DescriptorsTable descriptors = new DescriptorsTable(dtDescriptors, configuration, 0); 

		setStatus("Descriptors: " + (descriptors.getColumnSize() - descriptors.getConditionsSize()) + " conditions: " + descriptors.getConditionsSize());

		boolean supportMultilearning = configuration.getTrainingConfiguration() instanceof MultiLearningAbstractConfiguration;

		if(configuration.areMultiLearningData() && (!supportMultilearning || configuration.isForcedSingleTasklearning())){
			return  applyMultiModel(dtDescriptors, configuration);
		}

		DataTable dtResult = applyModel(descriptors, configuration);

		if (configuration.scaleY != null) // now we need to rescale results, if there is a scaling we use it
			configuration.scaleY.get().inverseScale(dtResult);

		return dtResult;
	}

	final private DataTable calculatePreFillTrain(final DataTable dtDescriptors, final DataTable dtValues, ModelAbstractConfiguration config) throws Exception{

		CompressedObject<Object> featureNetConfig = null;

		if(config instanceof ValidationConfiguration) {


			if(config.isFeatureNet()) 
			{
				// first step is to use STL and predict missed values
				ValidationConfiguration  stlConfig = (ValidationConfiguration) config.getDeepCopy();
				stlConfig.taskConfiguration.setFeatureNet(false);
				stlConfig.setFeatureNet(false);
				if(stlConfig.taskConfiguration instanceof MultiLearningAbstractConfiguration)
					((MultiLearningAbstractConfiguration)stlConfig.taskConfiguration).noMultiLearning = true;

				DataTable dtInterimResults = calculateTrain(dtDescriptors, dtValues, stlConfig); // will calculate STL predictions

				dtInterimResults.keepResults();

				if(dtInterimResults.getColumnsSize() != config.getOptions().size())
					throw new CriticalException("For feature the number of predictions " + dtInterimResults.getColumnsSize() + " != properties " + config.getOptions().size());

				dtDescriptors.mergeColumnsWith(dtInterimResults,QSPRConstants.FEATURE_NET_PREFIX);
				featureNetConfig = new CompressedObject<Object>(stlConfig.getOptimisedConfiguration());
			}

			else

				if(config.getTrainingConfiguration() instanceof LSSVMGConfiguration && config.areMultiLearningData() 
						&& !config.isForcedSingleTasklearning() && !(config instanceof BaggingConfiguration)) {

					final int trainingSetSize= dtValues.getRowsSize();
					final DescriptorsTable trainingSetDescriptors = new DescriptorsTable(dtDescriptors, config, trainingSetSize); // only training set molecules
					final LabelsTable values = new LabelsTable(dtValues, config); // duplicate

					Map<String,Integer[]> moleculesPerHash  = trainingSetDescriptors.groupMoleculesByMD5Hashes(values);
					boolean valuesArePresentForAllMolecules = true;
					for(String hash: moleculesPerHash.keySet()) {
						Integer vals[] = moleculesPerHash.get(hash);
						int n = vals[0]; 
						for(int m: vals)
							if(m != n) valuesArePresentForAllMolecules = false;
					}

					if(!valuesArePresentForAllMolecules) {
						// first step is to use STL and predict missed values

						DataTable extendedData = dtDescriptors.getSlice(0, trainingSetSize);
						ValidationConfiguration stlConf = (ValidationConfiguration) config.getDeepCopy();

						((MultiLearningAbstractConfiguration)stlConf.getTrainingConfiguration()).noMultiLearning = true; //STL
						DataTable dtInterimResults = calculateTrain(extendedData, dtValues, stlConf); // will calculate STL predictions, CV results and predictions will be used to feed models

						int steps = 3;
						for(int a = 1; a <= steps; a++) {
							System.out.println(dtInterimResults.getRowsNoErrorsSize()+" "+ dtInterimResults.getRowsSize());
							//next steps are to use predicted values to further improve predictions
							extendedData = dtDescriptors.getSlice(0, trainingSetSize); // forming new descriptor dataset
							DataTable extendedValues = dtValues.getSlice(0, trainingSetSize); // forming new labels dataset
							moleculesPerHash  = trainingSetDescriptors.groupMoleculesByMD5Hashes(values);
							for(int i=0;i<trainingSetSize;i++) { // adding additional rows
								String md5= trainingSetDescriptors.getMoleculeMD5(i, true);
								Integer vals[] = moleculesPerHash.get(md5);
								int max = Collections.max(Arrays.asList(vals));
								int start = extendedValues.getRowsSize();
								boolean error = false;
								for(int k=0;k<vals.length;k++) { // adding new data and values
									for(;vals[k]<max;vals[k]++){
										AbstractDataRow r = extendedData.getRow(i);
										extendedData.addRow(r.getDeeperCopy()); // adding a new datapoint to have full matrix
										extendedValues.addRow(); // adding a new value to have full matrix
										Serializable val = dtInterimResults.getValue(i, QSPRConstants.PREDICTION_RESULT_COLUMN+k);
										extendedValues.setValue(val);
										if(val==null || dtInterimResults.getRow(i).isError()) {
											error = true;
											break;
										}
										extendedValues.setValue(QSPRConstants.CLASS, k);
									}
								}

								if(error) // setting all records to error
									for(extendedValues.getRow(i).setError("calculations of submodel failed");start<extendedValues.getRowsSize();start++)
										extendedValues.getRow(start).setError("calculations of submodel failed");
							}

							extendedValues.removeErrors(); // clean from erroneous rows, if any

							if(a == steps) { //adding the test set for the last run
								int interimTrainingSetSize= extendedValues.getRowsSize();
								for(int i = trainingSetSize; i<dtDescriptors.getRowsSize();i++)
									extendedData.addRow(dtDescriptors.getRow(i));
								dtInterimResults = calculateTrain(extendedData, extendedValues, config); // will calculate single prediction using all data
								DataTable finalRes = dtInterimResults.getSlice(0, trainingSetSize);  // training part
								for(int i = interimTrainingSetSize; i<dtInterimResults.getRowsSize();i++) // restoring the part with predictions for the test set
									finalRes.addRow(dtInterimResults.getRow(i));
								return finalRes;
							}else 
								dtInterimResults = calculateTrain(extendedData, extendedValues, config.getDeepCopy()); // will calculate single prediction for each

							dtInterimResults.keepResults();
						}
					}
				}



		}

		DataTable res =  calculateTrain(dtDescriptors, dtValues, config);
		if(featureNetConfig != null)
			((ValidationConfiguration)config).taskConfiguration.featureNetConfig = featureNetConfig;

		return res;
	}

	void checkMinimaRecords(ModelAbstractConfiguration configuration, int records) {

		if(isRunningTest())return;
		int factor = configuration instanceof BaggingConfiguration?4:2;
		if(records<configuration.requireMinimumRecords()/factor)throw new 
		CriticalException("Implementation of this method cannot work with less than N="+
				configuration.requireMinimumRecords()/factor+" samples (this is a minimal number just used for debugging purposes) while the current set has only n="+records+" samples. It is impossible to calculate a model with so few data.");
	}

	final private DataTable calculateTrain(DataTable dtDescriptors, DataTable dtValues, ModelAbstractConfiguration configuration) throws Exception
	{
		int trainingSetSize =  dtValues.getRowsSize();
		if (trainingSetSize == 0)
			throw new CriticalException("Calculation of descriptors for all molecules failed. Calculation is not possible.");

		int num = configuration.requireMinimumRecords()/4;

		if(!configuration.skipSizeCheck() && (num  > trainingSetSize && !isRunningTest()))
			throw new CriticalException("This method " + configuration.getInformativeName() + 
					" require at least " + num + " records (so few records are used for debugging only), but after descriptor calculation only " + trainingSetSize + " were provided. Calculation is not possible.");

		/**
		 * All rows will contain information whether they are training or have validation type 
		 */
		for (int i = 0; i < dtDescriptors.getRowsSize(); i++)
		{
			AbstractDataRow r = dtDescriptors.getRow(i);
			r.removeAttachment(QSPRConstants.VALIDATION);
			if (i >= trainingSetSize)
				r.addAttachment(QSPRConstants.VALIDATION, true);
		}


		boolean supportMultilearning = configuration instanceof MultiLearningAbstractConfiguration || 
				(configuration instanceof ValidationConfiguration && ((ValidationConfiguration)configuration).taskConfiguration instanceof MultiLearningAbstractConfiguration);

		if(configuration.areMultiLearningData() && (!supportMultilearning || configuration.isForcedSingleTasklearning()))
			return  trainMultiModel(dtDescriptors, dtValues, configuration);  // multiple single models

		DescriptorsTable descriptors = new DescriptorsTable(dtDescriptors, configuration, trainingSetSize); 

		setStatus("Descriptors: " + (descriptors.getColumnSize() - descriptors.getConditionsSize()) + " conditions: " + descriptors.getConditionsSize());

		// create scaling and fill in it in configuration 

		DataTable dtResult = null;

		iterations = null;

		// initial values; scaling is added if it is in cfg
		LabelsTable values = new LabelsTable(dtValues, configuration);

		if(configuration.getExtendedDataType() == ModelAbstractConfiguration.DataType.REGRESSION && !configuration.isSupportRegression())
			throw new CriticalException("The selected method: " + supportedTaskType + " does not support regression tasks.");

		if (!configuration.isUploadedConfiguration()){
			checkMinimaRecords(configuration, values.getDataSize());
			dtResult = trainModel(descriptors, values, configuration);
			if (configuration.saveModels() && configuration instanceof DataDrivenConfiguration)
				((DataDrivenConfiguration) configuration).setTrainingSet(getTrainingSet(dtDescriptors, dtValues));
		}
		else
			dtResult = processUploadedModel(descriptors, values, configuration);

		if (configuration.scaleY != null) // now we need to rescale results, if there is a scaling we use it
			configuration.scaleY.get().inverseScale(dtResult);

		if(iterations != null) {
			configuration.iterations = iterations;
			setIterations(configuration);
		}
		return dtResult;

	}


	Object[] getReducedSets(DescriptorsTable dtDescriptors, DataTable dtExpValues, 
			ModelAbstractConfiguration configuration, int output) throws IOException{
		DataTable newVals = dtExpValues.getSlice(0, 0);
		DataTable newDesc = dtDescriptors.getRawData().getSlice(0, 0);

		Map<String,Integer> md5s = new HashMap<String,Integer>();

		for(int row = 0; row < dtExpValues.getRowsSize(); row++) {
			int outputnumber = (int) Math.rint((Double) dtExpValues.getValue(row, QSPRConstants.CLASS));

			if(output == outputnumber) { // only not errors should go there
				newVals.addRow(dtExpValues.getRow(row));
				newDesc.addRow(dtDescriptors.getRawRow(row));
				md5s.put(dtDescriptors.getMoleculeMD5(row, false),newDesc.currentRow); // all molecules with data are stored
			}
		}

		int missed = 0;

		if(configuration instanceof MultiLearningAbstractConfiguration) {
			MultiLearningAbstractConfiguration conf = (MultiLearningAbstractConfiguration) configuration;
			if(conf.getImplicit() != null && conf.getImplicit()[output] != null)
				for(int row = 0; row < dtExpValues.getRowsSize(); row++) {

					String md5 = dtDescriptors.getMoleculeMD5(row, false);

					if(!md5s.containsKey(md5)) { // really new molecule

						int outputnumber = (int) Math.rint((Double) dtExpValues.getValue(row, QSPRConstants.CLASS));
						if(output == outputnumber) throw new IOException("New values are found for old output"); // already added above

						missed++;
						AbstractDataRow r = dtExpValues.createRow();
						r.setValue(0, conf.getImplicit()[output]);
						newVals.addRow(r);
						newDesc.addRow(dtDescriptors.getRawRow(row));
						md5s.put(md5,newDesc.currentRow); // all molecules with data are stored
					}
				}

		}
		// adding all other molecules, which will be just predicted
		for(int row = 0; row < dtDescriptors.getDataSize(); row++) {

			String md5 = dtDescriptors.getMoleculeMD5(row, false);

			if(!md5s.containsKey(md5)) { // new molecule
				newDesc.addRow(dtDescriptors.getRawRow(row));
				md5s.put(md5,newDesc.currentRow); // all molecules with data are added
			}
		}

		setStatus("Created set with values: " + newVals.getRowsSize() + (missed >0? " including missed: " + missed : "") + " for property #" + output);

		ModelAbstractConfiguration config = configuration.getDeepCopy();
		config.scaleY = null; // we should re-init scaling using new data
		config.scaleX = null; 
		config.setTheOptions(new ArrayList<Integer>(Arrays.asList(new Integer[] { config.getOptions().get(output) })));

		Object tabs[] = new Object[4];

		tabs[0] = new DescriptorsTable(newDesc, config, newDesc.getRowsSize()); // create scaling and fill in it
		tabs[1] = new LabelsTable(newVals, config);
		tabs[2] = config;
		tabs[3] = md5s;

		return tabs;
	}

	private void trainOneProperty(int property,	DescriptorsTable dtDescriptors, DataTable dtExpValues,  ModelAbstractConfiguration configuration, ModelAbstractConfiguration configs[], DataTable predictions[]) throws Exception {
		Object sets[] = getReducedSets(dtDescriptors, dtExpValues, configuration, property);
		DescriptorsTable des = (DescriptorsTable) sets[0];
		LabelsTable lab = (LabelsTable) sets[1];
		ModelAbstractConfiguration config = (ModelAbstractConfiguration) sets[2];

		checkMinimaRecords(configuration, lab.getDataSize());
		DataTable res = trainModel(des, lab, config); // correctly validated results but not for all molecules ...
		predictions[property] = res.getSlice(0, 0);

		@SuppressWarnings("unchecked")
		Map<String,Integer> map = (Map<String,Integer>) sets[3];

		for(int row = 0; row < dtDescriptors.getDataSize(); row++) {
			String md5 = dtDescriptors.getMoleculeMD5(row, false);
			if(map.get(md5) == null)throw new IOException("Missed md5 " + md5);
			predictions[property].addRow(res.getRow(map.get(md5)));
		}
		configs[property] = config.getOptimisedConfiguration();

		if (!configuration.saveModels())
			configs[property].cleanBulkyStuff(); // do not store models if not required
		else
			if (configuration instanceof DataDrivenConfiguration)
				((DataDrivenConfiguration) configs[property]).setTrainingSet(getTrainingSet(des.getRawData().getSlice(0, lab.getDataSize()), lab.getRawData()));
	}

	DescriptorsTable getDescriptorsForFeatureNet(DataTable dtDescriptors,ModelAbstractConfiguration configuration,int property,int size) throws IOException {
		DataTable desc = dtDescriptors.getDeepCopy();
		HashSet <String> deleteColumns = new HashSet<String>();
		deleteColumns.add(desc.getColumn(dtDescriptors.getColumnsSize()-configuration.getOptions().size()+property)); // delete the property to be analyzed
		desc.deleteByList(deleteColumns);
		return new DescriptorsTable(desc, configuration, size); 
	}

	private DataTable trainMultiModel(DataTable dtDescriptors, DataTable dtExpValues, ModelAbstractConfiguration configuration) throws Exception{

		ModelAbstractConfiguration configs[] = new ModelAbstractConfiguration[configuration.getOptions().size()];

		DataTable predictions[] = new DataTable[configs.length];

		//DescriptorsTable dtDesc = dtDescriptors;

		DescriptorsTable descriptors = null;

		for(int property = 0; property < configs.length; property++) { // making model for each output and after that aggregating them

			if(configuration.isFeatureNet())
				descriptors=getDescriptorsForFeatureNet(dtDescriptors,configuration,property,dtExpValues.getRowsSize());
			else 
				if (descriptors == null) // first time only
					descriptors = new DescriptorsTable(dtDescriptors, configuration, dtExpValues.getRowsSize()); 

			setStatus("Descriptors: " + (descriptors.getColumnSize() - descriptors.getConditionsSize()) + " conditions: " + descriptors.getConditionsSize() + " property " + property + " ouf of " + configs.length);

			trainOneProperty(property, descriptors, dtExpValues, configuration, configs, predictions);
		}

		configuration = configuration.getOptimisedConfiguration();
		configuration.setModel(configs);

		return prepareCombinedTable(predictions);
	}

	private DataTable prepareCombinedTable(DataTable allres[]) {
		DataTable finalRes= new DataTable(true);

		boolean hasDM = allres[0].getColumnsSize() > 1 && allres[0].getColumn(1).contains(QSPRConstants.DM);

		for(int col = 0; col < allres.length ; col++) {
			finalRes.addColumn(QSPRConstants.PREDICTION_RESULT_COLUMN + col);
			if(hasDM) {
				String name = allres[0].getColumn(1).substring(QSPRConstants.DM.length());
				finalRes.addColumn(QSPRConstants.DM + col + name);
			}
		}

		for(int i = 0; i < allres[0].getRowsSize();i ++) { // developing a model for each output and after that aggregating them
			finalRes.addRow();
			int k =0;
			for(int j = 0; j < allres.length; j++) 
				if(allres[j].getRow(i).isError()) {
					finalRes.getRow(i).setError(allres[j].getRow(i).detailedStatus);
					break;
				}
				else {
					finalRes.setValue(i, k++, allres[j].getValue(i, 0));
					if(hasDM)
						finalRes.setValue(i, k++, allres[j].getValue(i, 1));
				}
		}

		return finalRes;
	}

	private DataTable applyMultiModel(DataTable dtDescriptors, ModelAbstractConfiguration configuration) throws Exception{

		if(configuration.isFeatureNet()) {
			DataTable dtInterimResults = applyMultiModel(dtDescriptors,(ModelAbstractConfiguration)configuration.featureNetConfig.get());
			dtInterimResults.keepResults();
			dtDescriptors.mergeColumnsWith(dtInterimResults,QSPRConstants.FEATURE_NET_PREFIX); // new set of descriptors
		}

		ModelAbstractConfiguration configs[] =  (ModelAbstractConfiguration[])configuration.getSavedModelAsObject();

		DataTable allres[] = new DataTable[configs.length];

		for(int property =0; property < configs.length; property++) {
			setStatus("applying property " + property + " out of " + configs.length);
			DescriptorsTable desctab = configuration.isFeatureNet()?getDescriptorsForFeatureNet(dtDescriptors, configuration, property, dtDescriptors.getRowsSize()):new DescriptorsTable(dtDescriptors, configs[property], 0); // create scaling and fill in it
			allres[property] = applyModel(desctab, configs[property]);
			if (configs[property].scaleY != null) // now we need to rescale results, if there is a scaling we use it
				configs[property].scaleY.get().inverseScale(allres[property]);
		}

		return prepareCombinedTable(allres);
	}


	/**
	 * Add method specific functions to init calculations
	 * 
	 * @throws InterruptedException
	 * @throws IOException
	 */
	void init() throws IOException, InterruptedException
	{

	}

	/**
	 * Cleaning (elimination of directories, making memory free) after the calculations
	 * 
	 * @throws InterruptedException
	 * @throws IOException
	 */

	void cleanup() throws IOException, InterruptedException
	{

	}

	private WorkflowNodeData getTrainingSet(DataTable dtDescriptors, DataTable dtExpValues) throws IOException
	{
		// Get the size of training set
		int nSample = dtExpValues.getRowsSize();
		WorkflowNodeData wnTrainingSet = new WorkflowNodeData();
		wnTrainingSet.addPort(dtDescriptors.getSlice(0, nSample));
		wnTrainingSet.addPort(dtExpValues.getSlice(0, nSample));
		return wnTrainingSet;
	}

	// The function is under construction, in testing! / Midnighter on Apr 12, 2011
	// I guess this function require re-implementation // IVT 30.08.2012
	// It is unlikely that it will work

	private void calculateDMs(WorkflowNodeData wndInput, DataTable dtPredictions, ModelAbstractConfiguration conf, boolean training) throws Exception
	{

		if (conf.dmModels == null)
			return;

		Set<String> dmNames = conf.dmModels.keySet();
		if (dmNames.isEmpty())
			return;

		DataTable dtDescriptors = wndInput.ports.get(0);

		// check whether structures were provided for AD
		DataTable dtStructures = null;

		int structuresPort = training ? 2 : 1;
		if (wndInput.ports.size() > structuresPort)
		{
			dtStructures = wndInput.ports.get(structuresPort);
			out.println("Structures were provided for AD calcualtions");
		}

		setStatus("Preparing to calculate " + dmNames.size() + " DMs: " + dmNames);
		DataTable dtSelectedDescriptors = getSelectedDescriptors(dtDescriptors, conf.getSelectedDescriptors());

		// Prepare a set of tasks
		CalculationTaskSet taskSet = new CalculationTaskSet(this);
		taskSet.prefix = "AD calculation: ";

		for (String dmName : dmNames)
		{
			CalculationTask cTask = new CalculationTask();
			int separator = dmName.indexOf("#");
			if (separator < 0)
				separator = dmName.length();
			cTask.taskName = dmName.substring(0, separator); // first part before # should contain ServerName to calculate DM
			cTask.configuration = conf.dmModels.get(dmName) != null ? conf.dmModels.get(dmName).get() : null;
			cTask.wndInput = new WorkflowNodeData(dtSelectedDescriptors);
			if (dtStructures != null)
				cTask.wndInput.addPort(dtStructures);

			taskSet.addTask(cTask);
		}

		// Run the calculation
		setStatus("Calculating DMs");
		taskSet.post();
		taskSet.calculate(false);

		// Process the results
		setStatus("Processing DM results");
		for (CalculationTask cTask : taskSet.tasks)
		{
			DataTable dtDmValues = cTask.getWndOutput().ports.get(0);
			String dmName = dtDmValues.getColumn(0);

			if (training) // store DM configurations for further application
			{
				Serializable dmModel = cTask.getWndOutput().ports.get(1).getValue(0, 0);
				conf.addDM(cTask.taskName, dmModel);
			}

			dtPredictions.addColumn(QSPRConstants.DM_PREFIX + dmName);
			dtDmValues.reset();
			while (dtDmValues.nextRow())
				dtPredictions.setValue(dtDmValues.currentRow, QSPRConstants.DM_PREFIX + dmName, dtDmValues.getValue());
		}
	}

	private void calculateDMsSafely(WorkflowNodeData wndInput, DataTable dtPredictions, ModelAbstractConfiguration conf, boolean training)
	{
		try
		{
			calculateDMs(wndInput, dtPredictions, conf, training);
		} catch (Exception e)
		{
			// We do not want to fail the whole model only because the DM calculation failed
			out.println("WARNING: could not calculate DMs. An exception occured: " + OCHEMUtils.exceptionToString(e));
		}
	}

	private DataTable getSelectedDescriptors(DataTable dtAllDescriptors, List<Integer> selectedDescriptorNums)
	{
		if (selectedDescriptorNums == null)
			return dtAllDescriptors;

		DataTable dtSelectedDescriptors = new DataTable(true);
		for (Integer descNum : selectedDescriptorNums)
			dtSelectedDescriptors.addColumn(dtAllDescriptors.getColumn(descNum));

		for (int i = 0; i < dtAllDescriptors.getRowsSize(); i++)
		{
			AbstractDataRow r = dtSelectedDescriptors.addRow();
			int col = 0;
			for (Integer descNum : selectedDescriptorNums)
				r.setValue(col++, dtAllDescriptors.getValue(i, descNum));
			if (dtAllDescriptors.getRow(i).getAttachment(QSPRConstants.VALIDATION) != null)
				r.addAttachment(QSPRConstants.VALIDATION, true);
		}
		return dtSelectedDescriptors;

	}

	void postTask(CalculationTaskSet taskSet, CalculationTask bagTask, boolean locally) throws IOException, ClassNotFoundException, InterruptedException
	{
		bagTask.lazyTaskRetrieval = true;
		taskSet.tasks.add(bagTask);
		bagTask.post(locally); // Specify, whether task will be calculated on the local server
		bagTask.wndInput = null;
	}

	public MachineLearningAbstractServer()
	{
		// Default flows for any machine learning servers
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setInputFlowGroup(2);
		setOutputFlowGroup(0);
	}


}
