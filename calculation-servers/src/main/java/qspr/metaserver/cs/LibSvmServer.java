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

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import com.eadmet.exceptions.CriticalException;

import qspr.metaserver.configurations.LibSvmConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LIBSVMUtils;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.util.AbstractSvmModel;
import qspr.metaserver.util.LibSvmModel;
import qspr.metaserver.util.SvmModel;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.CalculationTask;
import qspr.workflow.utils.CalculationTaskSet;
import qspr.workflow.utils.QSPRConstants;

public class LibSvmServer extends IgorsAbstractServer
{
	protected static String trainExe = "";
	protected static String predictExe = "";

	private final static String DATAFILE = "data_train";
	private final static String APPLYFILE = "data_test";
	private final static String MODEL = "model";

	private final static String gsTrainFileName = "data_gs_train";
	private final static String gsTestFileName = "data_gs_test";
	private final static String gsModelFileName = "model_gs";
	private final static String gsPredFileName = "data_gs_predictions";

	private final static String predictionsFileName = "data_predictions";

	private int realGridSteps = 0;

	public final static int MINUMIM_GRID_SIZE = 100;

	@Override
	public void setParam(String name, String value)
	{
		super.setParam(name, value);

		if ("EXE".equals(name))
			trainExe = value;
		if ("PREDICT_EXE".equals(name))
			predictExe = value;

	}

	private String deAlias(String value) {
		File f = new File(value);
		return getAliasedPath()+"/" + f.getName();
	}

	private String[] getLibSvmCommandLineTrain(LibSvmConfiguration conf, String file_data, String file_model)
	{
		List<String> list = new ArrayList<String>();

		list.add(deAlias(trainExe));
		list.add("-s");
		list.add(conf.svm_type + "");
		list.add("-t");
		list.add(conf.kernel_type + "");
		list.add("-d");
		list.add(conf.degree + "");
		if (conf.gamma != null)
		{
			list.add("-g");
			list.add(conf.gamma + "");
		}
		list.add("-r");
		list.add(conf.coef0 + "");
		if (conf.cost != null)
		{
			list.add("-c");
			list.add(conf.cost + "");
		}
		list.add("-n");
		list.add(conf.nu + "");
		if (conf.svrEpsilon != null)
		{
			list.add("-p");
			list.add(conf.svrEpsilon + "");
		}
		list.add("-m");		
		list.add(""+(getMemoryForExecutable()<100?100: getMemoryForExecutable())); // LibSVM memory is set to be half of all available to the server
		list.add("-e");
		list.add(conf.epsilon + "");
		list.add("-h");
		list.add("1");
		list.add("-b");
		list.add("0");

		if (conf.areClassificationData() && (conf.useWeighting != null && conf.useWeighting))
		{
			list.add("-w0");
			list.add("" + Math.round(10.0 * conf.classWeightRatio));
			list.add("-w1");
			list.add("10");
		}

		//		list.add("-q");
		list.add(getAliasedFileName(file_data));
		list.add(getAliasedFileName(file_model));

		return list.toArray(new String[0]);
	}

	private String[] getLibSvmCommandLineTest(LibSvmConfiguration conf, String file_data, String file_model, String file_predictions)
	{
		String[] commands = { deAlias(predictExe), file_data, file_model, file_predictions };
		return commands;

	}

	private void readPredictions(DataTable results, String file_predictions) throws Exception
	{
		BufferedReader br = getAliasedBufferedReader(file_predictions);
		String s = null;
		while ((s = br.readLine()) != null)
		{
			results.addRow();
			results.setValue(new Double(s));
		}
		br.close();
	}

	private double getBalancedAccuracy(LabelsTable real, DataTable predicted, String oneClassLabel) throws IOException
	{
		int classes = real.getNumberOfClasserPerProperty(0);
		int[][] confusion = new int[classes][classes];
		int[] count = new int[classes];
		real.getRawData().reset();
		predicted.reset();
		while (real.getRawData().nextRow())
		{
			predicted.nextRow();
			int realClass = Long.valueOf(Math.round((Double)real.getRawData().getValue())).intValue();

			int predClass;

			if (oneClassLabel != null) //We are one-class, +1 and -1 and stuff
			{
				if  (Long.valueOf(Math.round((Double)predicted.getValue())).intValue() == 1)
					predClass = Integer.valueOf(oneClassLabel);
				else
					predClass = (Integer.valueOf(oneClassLabel) + 1) % classes;
			}
			else
				predClass = Long.valueOf(Math.round((Double)predicted.getValue())).intValue();

			confusion[realClass][predClass]++;
			count[realClass]++;
		}
		double res = 0;
		for (int i=0; i<classes; i++)
			res += confusion[i][i] * 1.0 / count[i];
		res /= classes;
		return res;
	}

	private double getRMSE(LabelsTable real, DataTable predicted) throws IOException
	{
		real.getRawData().reset();
		predicted.reset();
		double sum = 0;
		int n = 0;
		while (real.getRawData().nextRow())
		{
			predicted.nextRow();
			double realValue = (Double)real.getRawData().getValue();
			double predValue = (Double)predicted.getValue();
			sum += (realValue - predValue) * (realValue - predValue);
			n++;
		}
		return Math.sqrt(sum / n);
	}




	protected double getCrossvalidatedAccuracy(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, LibSvmConfiguration receivedConf, int folds) throws Exception
	{
		String[] commands;
		double[] results = new double[folds];
		for (int fold = 0; fold < folds; fold++)
		{
			DescriptorsTable gsTrainingDescriptors = (DescriptorsTable)dtDescriptors.getSlice(0, 0);
			LabelsTable gsTrainingLabels = (LabelsTable)dtExpValues.getSlice(0, 0);
			DescriptorsTable gsTestDescriptors = (DescriptorsTable)dtDescriptors.getSlice(0, 0);
			LabelsTable gsTestLabels = (LabelsTable)dtExpValues.getSlice(0, 0);

			int foldSize = dtDescriptors.getDataSize() / folds;
			for (int i=0; i<dtDescriptors.getDataSize(); i++)
			{
				if (i >= fold * foldSize && i < (fold + 1) * foldSize)
				{
					gsTestDescriptors.getRawData().addRow(dtDescriptors.getRawRow(i));
					gsTestLabels.getRawData().addRow(dtExpValues.getRawRow(i));
				} else
				{
					gsTrainingDescriptors.getRawData().addRow(dtDescriptors.getRawRow(i));
					gsTrainingLabels.getRawData().addRow(dtExpValues.getRawRow(i));
				}
			}

			String gsTrainFold = gsTrainFileName + "_" + fold; 
			String gsTestFold = gsTestFileName + "_" + fold; 
			
			LIBSVMUtils.writeLibSvmTrainingSet(gsTrainingDescriptors, gsTrainingLabels, receivedConf.oneClass ? receivedConf.oneClassLabel : null,  getAliasedFileName(gsTrainFold));
			commands = getLibSvmCommandLineTrain(receivedConf, gsTrainFold, gsModelFileName);
			executeBinary(commands, stdout);

			LIBSVMUtils.writeLibSvmTestSet(gsTestDescriptors, getAliasedFileName(gsTestFold));
			commands = getLibSvmCommandLineTest(receivedConf, gsTestFold, gsModelFileName, gsPredFileName);
			executeBinary(commands, gsPredFileName);

			DataTable predictions = new DataTable(true);
			predictions.addColumn("RESULT");
			readPredictions(predictions, gsPredFileName);

			(new File(getAliasedFileName(gsModelFileName))).delete();
			(new File(getAliasedFileName(gsPredFileName))).delete();

			if (receivedConf.areClassificationData())
				results[fold] = getBalancedAccuracy(gsTestLabels, predictions, receivedConf.oneClass ? receivedConf.oneClassLabel : null);
			else
				results[fold] = getRMSE(gsTestLabels, predictions);
		}


		double res = 0;
		if (receivedConf.areClassificationData())
		{
			for (int i=0; i < folds; i++)
				res += results[i];
			return res / folds;
		} else
		{
			for (int i=0; i < folds; i++)
				res += results[i] * results[i];
			return -Math.sqrt(res / folds); // "minus" so that bigger == better
		}
	}

	protected LibSvmConfiguration gridSearch(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, LibSvmConfiguration receivedConf) throws Exception
	{
		//Form GridSearch subset
		setStatus("Performing grid search parameter optimization");

		DescriptorsTable desc = dtDescriptors;
		LabelsTable labels = dtExpValues;

		if(receivedConf.stepStart == null){
			int gridSearchSetSize = (int) Math.max(Math.round(receivedConf.gridSearchSetSize * dtExpValues.getDataSize()), Math.min(MINUMIM_GRID_SIZE, dtExpValues.getDataSize()));

			if (gridSearchSetSize <= 10) // to small? using all rows!
				throw new CriticalException("The size of grid search dataset too small (" + gridSearchSetSize + " rows). Use LibSVM with larger datasets.");

			desc = (DescriptorsTable) dtDescriptors.getSlice(0, 1);
			labels = (LabelsTable) dtExpValues.getSlice(0, 1);

			List<Integer> rows = new ArrayList<Integer>();

			for (int i=1; i<dtExpValues.getDataSize(); i++) 
				rows.add(i);

			Collections.shuffle(rows, new Random(receivedConf.getSeed()));
			rows = rows.subList(0, gridSearchSetSize - 1);

			for (Integer row : rows) 
			{
				desc.getRawData().addRow(dtDescriptors.getRawRow(row));
				labels.getRawData().addRow(dtExpValues.getRawRow(row));
			}	

			if(receivedConf.PARALLEL != null)
				return calculateInParallel(desc,labels,receivedConf);
		}

		realGridSteps = receivedConf.stepStart == null ? receivedConf.totalGridSteps() : receivedConf.stepStop - receivedConf.stepStart;

		receivedConf.bestAccuracy = Double.NEGATIVE_INFINITY;

		resetStartTime();

		int currentStep = 0;

		int folds = 3;
		
		for(Integer step: receivedConf.prepareSteps())
		{
			currentStep++;

			if(receivedConf.stepStart != null && currentStep < receivedConf.stepStart)continue;
			if(receivedConf.stepStop != null && currentStep == receivedConf.stepStop)break;

			receivedConf.configurationForGridStep(step);

			double currentAcc = getCrossvalidatedAccuracy(desc, labels, receivedConf, folds);

			if (currentAcc > receivedConf.bestAccuracy)
			{
				receivedConf.bestAccuracy = currentAcc;
				receivedConf.saveBestParameters();
				out.println("improved accuracy " + Math.abs(receivedConf.bestAccuracy) + "for "+receivedConf.getActualParametersString());
			}

			int elapsedStep = receivedConf.stepStart != null? currentStep - receivedConf.stepStart : currentStep;

			setPercentageCompleted(1.0 * elapsedStep / (realGridSteps + 1));
			setStatus("Performing grid search parameter optimization (step "+elapsedStep+" out of "+realGridSteps+"): " + 
					receivedConf.getActualParametersString() + " " + " to finish in "	+ getTimeToComplete());
		}

		receivedConf.restoreBestParameters();

		out.println("The BEST configuration according to the grid search:");
		out.println(receivedConf.getActualParametersString());
		out.println("The accuracy of the best config is " + Math.abs(receivedConf.bestAccuracy));

		
		for (int fold = 0; fold < folds; fold++)
		{
			(new File(getAliasedFileName(gsTrainFileName + "_" + fold))).delete();
			(new File(getAliasedFileName(gsTestFileName + "_" + fold))).delete();
		}

		return receivedConf;
	}


	private LibSvmConfiguration calculateInParallel(
			DescriptorsTable dtDescriptorsOld, LabelsTable dtExpValuesOld,
			LibSvmConfiguration receivedConf) throws Exception {

		setStatus("Starting parallel calculations ...");

		DescriptorsTable dtDescriptors=(DescriptorsTable)dtDescriptorsOld.getSlice(0, dtExpValuesOld.getDataSize()); // we do not need to send test set for calculations

		int totalSteps = receivedConf.totalGridSteps();

		int step= totalSteps / receivedConf.PARALLEL <1 ? 1 : totalSteps / receivedConf.PARALLEL;

		Boolean saveModel = receivedConf.saveModels; // to preserve the stage
		receivedConf.saveModels = false; // we do not save sub models in parallel calculations 

		CalculationTaskSet taskSet = new CalculationTaskSet(this);
		for(int n = 0 ; n < totalSteps ; n += step){

			receivedConf.stepStart = n;
			receivedConf.stepStop = n + step;

			CalculationTask cTask = new CalculationTask();
			cTask.taskName = supportedTaskType; 

			cTask.configuration = receivedConf;
			cTask.setParent(this);
			cTask.wndInput = new WorkflowNodeData(dtDescriptors.getRawData()).addPort(dtExpValuesOld.getRawData());
			taskSet.tasks.add(cTask);
			cTask.lazyTaskRetrieval = true;
			cTask.post();
		}

		taskSet.calculate(true);

		LibSvmConfiguration best = (LibSvmConfiguration)taskSet.tasks.get(0).getWndOutput().ports.get(1).getValue(0, 0);

		for(int i=1;i<taskSet.tasks.size();i++){
			LibSvmConfiguration modelConfiguration = (LibSvmConfiguration)taskSet.tasks.get(i).getWndOutput().ports.get(1).getValue(0, 0);

			out.println("best accuracy " + best.bestAccuracy + " new accuracy: " + 
					modelConfiguration.bestAccuracy + "for "+modelConfiguration.getActualParametersString());

			if(modelConfiguration.bestAccuracy <= best.bestAccuracy) continue;

			best = modelConfiguration;
			out.println("***** best accuracy was improved!");
		}

		best.stepStart=best.stepStop=null;
		receivedConf.saveModels = saveModel;

		return best;
	}

	@Override
	protected DataTable applyModelBatch(DescriptorsTable dtDescriptors, ModelAbstractConfiguration receivedConf) throws Exception
	{
		LibSvmConfiguration conf = (LibSvmConfiguration) receivedConf;

		LIBSVMUtils.writeLibSvmTestSet(dtDescriptors, getAliasedFileName(APPLYFILE));

		AbstractSvmModel model = (AbstractSvmModel) conf.getSavedModelAsObject();
		model.writeModelToFile(getAliasedFileName(MODEL));

		String[] commands = getLibSvmCommandLineTest(conf, APPLYFILE, MODEL, predictionsFileName);

		executeBinary(commands, predictionsFileName);

		DataTable dtResults = new DataTable(true);
		dtResults.addColumn(QSPRConstants.PREDICTION_RESULT_COLUMN);

		readPredictions(dtResults, predictionsFileName);

		return dtResults;
	}

	private void autoDetectConfiguration(LibSvmConfiguration conf)
	{
		if (conf.areClassificationData())
		{
			if (conf.oneClass)
				conf.svm_type = 2;
			else if (!conf.useNu)
				conf.svm_type = 0; //C-SVC
			else
				conf.svm_type = 1;
		}
		else 
		{
			if (!conf.useNu)
				conf.svm_type = 3; //epsilon-SVR;
			else
				conf.svm_type = 4; //nu-SVR;
		}
	}

	private void fillDefaultParameters(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, LibSvmConfiguration configuration)
	{
		if (configuration.gamma == null)
			configuration.gamma = 1D / dtDescriptors.getColumnSize();
	}

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws Exception
	{
		LibSvmConfiguration conf = (LibSvmConfiguration) receivedConf;

		autoDetectConfiguration(conf);

		LIBSVMUtils.writeLibSvmTrainingSet(dtDescriptors, dtExpValues, conf.oneClass ? conf.oneClassLabel : null, getAliasedFileName(DATAFILE));

		System.gc();

		setStatus("Available memory buffer " + getMemoryForExecutable() + "MB");

		if (conf.gridSearch)
		{
			LibSvmConfiguration newConfig = gridSearch(dtDescriptors, dtExpValues, conf);
			setPercentageCompleted( (conf.PARALLEL == null ? 1: conf.PARALLEL) * realGridSteps / (realGridSteps + 10.));
			setStatus("Training model with the optimized configuration - to finish in " + getTimeToComplete());

			conf.copy(newConfig);

			if(conf.stepStart!=null){ // this is a submodel
				DataTable dtResults = new DataTable(true);
				dtResults.addColumn(QSPRConstants.PREDICTION_RESULT_COLUMN);
				for(int i=0;i<dtDescriptors.getDataSize();i++){
					dtResults.addRow();
					dtResults.setValue(0.);
				}
				return dtResults;
			}

		} else
		{
			fillDefaultParameters(dtDescriptors, dtExpValues, conf);
			setStatus("Training the model with the given configuration");
		}

		String[] commands = getLibSvmCommandLineTrain(conf, DATAFILE, MODEL);
		executeBinary(commands, MODEL);

		AbstractSvmModel svmModel = null;

		try{ // first trying Sergey's format, which provides models of a smaller size but cannot handle > 
			svmModel = new LibSvmModel();
			svmModel.readModelFromFile(getAliasedFileName(MODEL));
			conf.setModel(svmModel);
		}catch(IOException e){ // failed?! Trying the one based on DataTable
			svmModel = new SvmModel();
			svmModel.readModelFromFile(getAliasedFileName(MODEL));
		}

		conf.setModel(svmModel);

		DataTable dtResults = applyModel(dtDescriptors, conf);
		setStatus("SVM finished");
		return dtResults;
	}

	public LibSvmServer()
	{
		supportedTaskType = QSPRConstants.LIBSVM;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}

}
