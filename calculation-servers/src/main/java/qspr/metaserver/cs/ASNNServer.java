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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.Writer;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import qspr.metaserver.configurations.ASNNConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.cs.util.ASNNEnsemble;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.util.DataScaling;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.CalculationTask;
import qspr.workflow.utils.CalculationTaskSet;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.utils.FileUtils;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.OSType;
/**
 * 
 * @author itetko
 * Implements abstract layer server for several methods
 * developed under common C/C++ interface
 */
public class ASNNServer extends MultiLearningAbstractServer
{
	public static final String PARALLEL="PARALLEL";
	public static final String GRAPHICS="graphics";
	public static final String FRESHMODEL="models";
	public static final double MAX_PROB_STD=10;

	protected String VALUESSEPARATOR = "\t";

	String DATAFILE = "data", CFG = "cfg", RESULTS = "results.txt", MODEL = "model", OUTPUT = "output";

	@Override
	public DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws Exception
	{
		DescriptorsTable dtOriginal = dtDescriptors;
		ASNNConfiguration configuration=(ASNNConfiguration)receivedConf;

		if(configuration.ensemble>256) throw new CriticalException("ASNN can work only with 256 or less neural networks in the ensemble");

		//store some critical values
		int iterations=configuration.iterations; // we need to preserve them
		boolean asnn=configuration.asnn; // we also preserve asnn
		if(configuration.libraryModel!=null)
			prepareLibrary(dtDescriptors, dtExpValues, configuration);
		else
		{
			// check whether parallelization is required
			if(configuration.additionalParam != null){
				Map<String,String> params=getAdditionalParam(configuration.additionalParam);
				if(params.containsKey(PARALLEL)){
					if(params.containsKey("PRUNE"))throw new IOException("PRUNE and PARALLEL are not compatible: "+configuration.additionalParam);
					int n= (int)Long.parseLong(params.get(PARALLEL));
					if(configuration.ensemble>n && n!=0){
						configuration.asnn=false; // no asnn for submodels  -- to speed-up analysis
						configuration=calculateInParallel(dtDescriptors, dtExpValues, configuration,Integer.parseInt(params.get(PARALLEL)));
						configuration.createASNNApplyModelConfiguration();
						configuration.asnn=asnn; // restore asnn
					}
				}
			}
		}

		saveAggregatedData(getAliasedFileName(DATAFILE),dtDescriptors, dtExpValues, receivedConf);

		createCfgFile(configuration,dtExpValues.getColumnSize());

		setStatus("Starting ASNN calculations ...");
		String[] commands = new String[] { getExeFile(), CFG, DATAFILE };
		executeBinary(commands, FRESHMODEL, 0);
		configuration=readConfiguration(configuration); // we need to get parameters of the resulting cfg
		configuration.iterations=iterations; // preserve iterations; their number was changed to 1 in case of parallel calculations

		ASNNEnsemble n=new ASNNEnsemble();
		BufferedReader model = getAliasedBufferedReader(FRESHMODEL);
		n.readModelNN(model);
		model.close();

		configuration.descriptors = loadDescriptors(dtDescriptors.columnNames);

		if(configuration.descriptors != null){ // we have to re-save only filtered descriptors on this step; this was pruning
			setStatus(""+configuration.descriptors.size()+" descriptors were selected: "+configuration.descriptors.toString());
			saveAggregatedData(getAliasedFileName(DATAFILE),dtDescriptors, dtExpValues, receivedConf);
		}

		if(configuration.asnn){
			setStatus("Preparing ASNN ranks ...");
			n.prepareRanksAndBiases(configuration,getAliasedFileName(DATAFILE)); // getting ranks for the training set but only for the ANN
		}
		//compress cfg, data and model and encode it for later use in applier mode

		configuration.setModel(FileUtils.getFileAsString(getAliasedFileName(FRESHMODEL))); // stored as String ... old incompatibility problem

		// now adding predictions for training and test sets
		return applyModel(dtOriginal, configuration); // now we apply model to all data (both training and test sets)

	}


	private LinkedHashMap<String,Integer> loadDescriptors(List<String>columns) throws IOException {
		int descSize = columns.size();
		BufferedReader descriptors = getAliasedBufferedReader("pruned");
		String s=descriptors.readLine();
		if(s.length()!=descSize)throw new IOException(" The number of descriptors in pruned "+
				s.length()+"!= number "+descSize+" of columns!");

		LinkedHashMap<String,Integer> columnsToKeep = new LinkedHashMap<String,Integer>();

		for(int i=0;i<s.length();i++){
			if(Integer.valueOf(""+s.charAt(i))==0){
				columnsToKeep.put(columns.get(i),i);
			}
		}

		if(columnsToKeep.size()==descSize)return null;

		return columnsToKeep;
	}


	@Override
	protected DataTable applyModelMulti(DescriptorsTable dtDescriptors, MultiLearningAbstractConfiguration receivedConf) throws IOException
	{
		ASNNConfiguration configuration=(ASNNConfiguration)receivedConf;

		ASNNEnsemble n=new ASNNEnsemble();
		BufferedReader model = new BufferedReader(new StringReader(configuration.getSavedModelAsObject().toString()));
		n.readModelNN(model);
		model.close();

		int outputs = n.outputs();

		saveAggregatedData(getAliasedFileName(DATAFILE),dtDescriptors, null, receivedConf);

		if(configuration.additionalParam != null && configuration.additionalParam.contains("MEDIAN"))
			throw new CriticalException("Very old config cannot be applied");

		n.predictASNN(getAliasedFileName(DATAFILE), getAliasedFileName(GRAPHICS), configuration, this);

		return  readResult(GRAPHICS, outputs, configuration.getOptions(),configuration.scaleY==null?1:2);
	}


	/**
	 *  Parallelize calculations of the ASNN task on different servers
	 * @param dtDescriptors
	 * @param dtExpValues
	 * @param configuration
	 * @param networksNumber
	 * @throws Exception 
	 */


	protected ASNNConfiguration calculateInParallel(final DescriptorsTable dtDescriptorsOld, final LabelsTable dtExpValuesOld, ASNNConfiguration configuration, int networksNumber) throws Exception{

		setStatus("Starting parallel calculations ...");

		int ensemble=configuration.ensemble;
		boolean asnn=configuration.asnn;
		int seed=configuration.getSeed();

		configuration.models=1;  // to save data
		configuration.asnn=false;    // to avoid ASNN correction

		DescriptorsTable dtDescriptors=(DescriptorsTable)dtDescriptorsOld.getSlice(0, dtExpValuesOld.getDataSize()); // we do not need to send test set for calculations

		CalculationTaskSet taskSet = new CalculationTaskSet(this);
		for(int n=0;n<ensemble;n+=networksNumber){
			configuration.ensemble=	networksNumber < (ensemble-n)?networksNumber : (ensemble-n);  // to have exactly the same number of models

			CalculationTask cTask = new CalculationTask();
			cTask.taskName = "ASNNP"; // will be calculated as ASNNP
			configuration.seed = n+1; // new seed for each network

			cTask.configuration = configuration;
			cTask.setParent(this);
			cTask.wndInput = new WorkflowNodeData(dtDescriptors.getRawData()).addPort(dtExpValuesOld.getRawData());
			taskSet.tasks.add(cTask);
			cTask.lazyTaskRetrieval = true;
			cTask.post();
		}

		taskSet.calculate(false);

		StringBuffer previousmodel=new StringBuffer();

		for(int i=0;i<taskSet.tasks.size();i++){
			ASNNConfiguration modelConfiguration = (ASNNConfiguration)taskSet.tasks.get(i).getWndOutput().ports.get(1).getValue(0, 0);
			saveModel(modelConfiguration);
			BufferedReader f=getAliasedBufferedReader(MODEL);
			String line;
			int models=0;
			for(int linenumber=0;(line=f.readLine())!=null;linenumber++){

				if(linenumber<3){
					if(i==0)previousmodel.append(line+"\n");
					continue;
				}

				if(linenumber==3){ // we need it to determine number of models and put correct number for ensemble
					String []vals=line.split("\\s+");
					models=(int)Long.parseLong(vals[0]);

					if(i==0){
						vals[0]=""+ensemble;
						for(String s: vals){
							previousmodel.append(s+"\t");
						}
						previousmodel.append("\n");
					}
					continue;
				}

				previousmodel.append(line+"\n");
				if(linenumber==((3*models)+3))break;
			}
			f.close();
		}

		BufferedWriter f=getAliasedBufferedWriter(MODEL);
		f.write(previousmodel.toString());
		f.close();
		configuration.asnn=asnn;
		configuration.ensemble=ensemble; // restore initial ensemble size
		configuration.seed=seed; // restore default seed
		return configuration;
	}

	private void prepareLibrary(DescriptorsTable dtDescriptors, LabelsTable dtExpValues,
			ASNNConfiguration configuration) throws IOException, InterruptedException {

		if(dtExpValues.getColumnSize()!=1)throw new IOException("LIBRARY mode is determined only for only one output ");
		if(configuration.libraryOutput==null){
			configuration.libraryOutput=0; // first output is used by default
		}

		configuration.savedmodel=configuration.libraryModel;
		saveModel(configuration);
		updateModel(configuration.libraryOutput);

		configuration.iterations=0; // we just need to predict the values
		configuration.models=2; // i.e., we need to read the model file
		configuration.additionalParam=null; // no additional params
		configuration.asnn=false; // no asnn at this step!
		//create data file for ANN Training

		dtExpValues.setScaling(null);

		saveAggregatedData(getAliasedFileName(DATAFILE), dtDescriptors, dtExpValues, configuration);

		createCfgFile(configuration,1);

		setStatus("Starting ann calculations ...");
		// we will need to read output from graphics file
		String[] commands = new String[] { getExeFile(), CFG, DATAFILE };
		executeBinary(commands,GRAPHICS,0);

		if(configuration.iterations == 0) { // only if we keep old model

			setStatus("ANN has finished. Reading results ...");
			DataTable values = new DataTable(true);
			DataTable results = new DataTable(true);

			BufferedReader graphics = getAliasedBufferedReader(GRAPHICS);

			String line;
			while ((line = graphics.readLine()) != null)
			{
				line = line.trim();
				if (!line.startsWith("mol_"))continue;
				values.addRow();
				results.addRow();
				String[] predictions = line.split("\\s+");
				values.setValue("exp", predictions[1]);
				results.setValue("pred", predictions[2]);
			}

			graphics.close();

			DataScaling scale=new DataScaling();
			scale.normByLinearRegression(values, results);
			dtExpValues.setScaling(scale);
			configuration.addScaleY(scale); 

			cleanup(); // now we delete the previous directory
			init(); // and recreate it again to keep clean
			saveModel(configuration);
			updateModel(configuration.libraryOutput);
		}

		configuration.libraryModel=null;
		configuration.setModel(null); // cleaning up old model

		configuration.createASNNApplyModelConfiguration();
	}



	private void updateModel(Integer output) throws IOException {

		if(output==null)throw new IOException("The target utput is not determined");

		PrintWriter newModel = null;
		BufferedReader oldModel = null;

		try{
			File fNew=new File(getAliasedFileName(FRESHMODEL));
			File fOld=new File(getAliasedFileName(MODEL));
			newModel = new PrintWriter(new FileWriter(fNew));
			oldModel = new BufferedReader(new FileReader(fOld));

			int outputNumber=0,ensemble=0,layers=0,neurons=0,bias=0;
			String line,values[];
			// reading the header
			for(int n=0;n<4;n++)
			{
				if((line = oldModel.readLine()) == null)
					throw new IOException("Cannot read header of model file, line=" +n);
				values = line.split("\\s+");
				if(n<2)newModel.println(line);
				if(n==2){// save correct normalization for the output value
					newModel.println("1\t1\t0"); // we do not care, since re-normalization will be performed in any case 
					//					newModel.println("1\t"+values[1+output*2]+"\t"+values[2+output*2]);
					outputNumber=Integer.parseInt(values[0]);
					if(outputNumber<=output)
						throw new IOException("The number of outputs in model "+outputNumber+
								" less than required output number "+output);
				}
				if(n==3){
					ensemble=Integer.parseInt(values[0]);
					bias=Integer.parseInt(values[2]);
					layers=Integer.parseInt(values[3]); // number of layers
					neurons=Integer.parseInt(values[2+layers]); // neurons on the layer before the last
					for(int i=0;i<values.length;i++)newModel.print(""+(i==(3+layers)?"1":values[i])+"\t");
					newModel.println();
				}
			}
			layers--;
			neurons+=bias==0?0:1;
			// reading the model
			for(int i=0;i<ensemble;i++)
				for(int j=0;j<layers;j++){
					if((line = oldModel.readLine()) == null)
						throw new IOException("Cannot read model file network="+i+" layers="+(j+2));
					if(j!=layers-1)newModel.println(line);
					else{
						values = line.split("\\s+");
						newModel.print(""+neurons);
						for(int k=neurons*output+1;k<(neurons*(output+1)+1);k++){
							newModel.print("\t"+values[k]);
						}
						newModel.println();
					}
				}
			newModel.println();

			newModel.close();
			oldModel.close();
			fOld.delete();
			fNew.renameTo(fOld);
		}catch(IOException e){
			throw new CriticalException("Failed in reading model file:"+e.getMessage());
		}finally{
			if(newModel!=null)newModel.close();
			if(oldModel!=null)oldModel.close();
		}
	}



	// here correct calculation of probability can be used
	// actually just 1/z is enough, since probability value is proportional to it
	private static double Phi(double z) {
		z=1./z;
		return z>MAX_PROB_STD?MAX_PROB_STD:z;
	}

	private DataTable readResult(String predictionFile, int outputs,List<Integer> optionsNumber,int libraryMode) throws IOException
	{
		BufferedReader graphics = getAliasedBufferedReader(predictionFile);

		DataTable dtResult = new DataTable(true);
		dtResult.id = QSPRConstants.PREDICTION_RESULT_COLUMN;

		String line;
		while ((line = graphics.readLine()) != null)
		{
			line = line.trim();
			if (!line.startsWith("mol_"))continue;
			dtResult.addRow();
			String[] predictions = line.split("\\s+");

			int offset=1+outputs; // a starting column with predictions
			int classNumber=0;

			for (int output = 0; output < optionsNumber.size(); output++,offset+=classNumber){
				String name=optionsNumber.size()==1?"":""+output;

				classNumber=optionsNumber.get(output);
				double std=0;
				if(classNumber>1){ //we have a classification task

					int predictedClass=LabelsTable.getPredictedClass(offset,classNumber,predictions);

					if(classNumber>2)
						dtResult.setValue(QSPRConstants.PREDICTION_RESULT_COLUMN + name, predictedClass); // class for multi-class classification
					else
						dtResult.setValue(QSPRConstants.PREDICTION_RESULT_COLUMN + name, LabelsTable.parseMissedValue(predictions[offset+1])); // just a value for compatibility with other code

					std=Double.parseDouble(predictions[offset+predictedClass+outputs]);
					double lag=classLag(offset,classNumber,predictions);
					dtResult.setValue(QSPRConstants.DM + name + QSPRConstants.PROBSTD,Phi(lag/std));
					dtResult.setValue(QSPRConstants.DM + name + QSPRConstants.STDEV,std);
					dtResult.setValue(QSPRConstants.DM + name + QSPRConstants.CLASSLAG,1-lag);
					dtResult.setValue(QSPRConstants.DM + name + QSPRConstants.CORREL,
							Double.parseDouble(predictions[offset+predictedClass+outputs*2]));

				}else{
					dtResult.setValue(QSPRConstants.PREDICTION_RESULT_COLUMN + name, LabelsTable.parseMissedValue(predictions[offset]));
					//					out.println("set value "+name+" "+parseAnnValue(predictions[offset]));
					switch(libraryMode){
					case 1: // for normal calculations
						dtResult.setValue(QSPRConstants.DM + name + QSPRConstants.STDEV,
								Double.parseDouble(predictions[offset+outputs]));
						dtResult.setValue(QSPRConstants.DM + name + QSPRConstants.CORREL,
								Double.parseDouble(predictions[offset+outputs*2]));
						break;
					case 2: // for library mode
						dtResult.setValue(QSPRConstants.DM + name + QSPRConstants.CORREL,
								Double.parseDouble(predictions[offset+outputs*2]));
						dtResult.setValue(QSPRConstants.DM + name + QSPRConstants.STDEV,
								Double.parseDouble(predictions[offset+outputs]));
						break;
					default: break; // just avoid 
					}
				}
			}
		}
		graphics.close();
		return dtResult;
	}


	private ASNNConfiguration readConfiguration(ASNNConfiguration configuration) throws IOException{

		configuration.models = 2;
		configuration.k = 0;
		configuration.parzen = 0.;

		if(!configuration.asnn) return configuration; // no parameters are required to read

		List<String> additional = OCHEMUtils.getFileAsStringList(getAliasedFileName(OUTPUT));

		for (String string : additional) {
			if(string.contains("EARLY_STOPPING_POINT S1 ASNN sigma=")){
				String[] arrayPM = string.split(" ");
				for (String pm : arrayPM) {
					if(pm.contains("KNN=")){
						pm = pm.replace("KNN=", "").trim();
						configuration.k=(int)Long.parseLong(pm);
						continue;
					}
					if(pm.contains("sigma=")){
						pm = pm.replace("sigma=", "").trim();
						configuration.parzen=Double.parseDouble(pm);
					}


				}
				break;	
			}
		}

		if(configuration.k == 0){ // no ASNN correction is required
			configuration.k = 0;
			configuration.asnn = false; 
			configuration.parzen = 0.;
		}

		return configuration;
	}

	private void createCfgFile(ASNNConfiguration configuration, int outputs) throws IOException{

		out.println("Writing cfg file");

		// standard parameters
		String parameters = "FILE=data"+OSType.endLine()+
				"NORM=1"+OSType.endLine()+"NAMES=1"+OSType.endLine()+"ALLPOINTS=0";

		parameters += OSType.endLine()+"PRINT=1";
		parameters += OSType.endLine()+"NONZERO=0"; // no filtering will be used
		parameters += OSType.endLine()+"CORRELATION=0"; // no de-correlation will be used

		parameters += OSType.endLine()+"SEED="+configuration.getSeed();
		parameters += OSType.endLine()+"ASNN="+(configuration.asnn?"1":"0");
		parameters += OSType.endLine()+"ENSEMBLE="+configuration.ensemble;
		parameters += OSType.endLine()+"NEURONS="+configuration.neurons;
		parameters += OSType.endLine()+"TRAINING="+configuration.training;
		parameters += OSType.endLine()+"OUTPUTS="+outputs;
		parameters += OSType.endLine()+"ITERATIONS="+(configuration.iterations>0?configuration.iterations:0);
		parameters += OSType.endLine()+"MODELS="+configuration.models;
		parameters += OSType.endLine()+"K="+configuration.k;
		parameters += OSType.endLine()+"PARZEN="+configuration.parzen;
		parameters += OSType.endLine()+"MISSED="+LabelsTable.MISSED_VALUES; // just to avoid possible problems

		if(configuration.additionalParam != null){
			Map<String,String> params=getAdditionalParam(configuration.additionalParam);
			for(String name: params.keySet())
				parameters += OSType.endLine()+name+"="+params.get(name);
		}

		parameters += OSType.endLine()+"STOP"+OSType.endLine()+"END"+OSType.endLine();

		FileUtils.saveStringToFile(parameters, getAliasedFileName(CFG));
	}

	/**
	 * Return difference between predicted (largest) and second largest value
	 * @param offset
	 * @param classNumber
	 * @param predictions
	 * @return
	 */
	static private double classLag(int offset,int classNumber,String [] predictions){

		int predictedClass=LabelsTable.getPredictedClass(offset,classNumber,predictions,-1);
		int secondLargestClass=LabelsTable.getPredictedClass(offset,classNumber,predictions,predictedClass);

		return 	LabelsTable.parseMissedValue(predictions[predictedClass+offset])-
				LabelsTable.parseMissedValue(predictions[secondLargestClass+offset]); // this value is always non-negative
	}

	private void saveModel(ASNNConfiguration modelConf) throws IOException {
		// decode the data

		if(modelConf.savedmodel ==null)throw new IOException("model is absent");
		FileUtils.saveStringToFile(modelConf.getSavedModelAsObject().toString(), getAliasedFileName(MODEL));

	}

	/**
	 * Return Map with additional parameters, which are submitted as PARAM=VALUE
	 * pairs
	 * 
	 * @param additionalParam
	 * @return
	 */
	private Map<String, String> getAdditionalParam(String additionalParam) {
		Map<String, String> params = new HashMap<String, String>();

		for (String param : additionalParam.split(",")) {
			String[] values = param.split("=");
			if (values == null || values.length != 2)
				continue;
			params.put(values[0], values[1]);
		}
		return params;
	}

	@Override
	public int saveAggregatedData(String filename, DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws IOException {
		BufferedWriter bw = getAliasedBufferedWriter(filename);
		bw.write("" + saveAggregatedDescriptorsAndValues(dtDescriptors, dtExpValues,null,receivedConf).size() +"\n");
		List <Integer> savedRecords = saveAggregatedDescriptorsAndValues(dtDescriptors, dtExpValues, bw,receivedConf); // how many we can compress
		bw.close();

		ASNNConfiguration conf = (ASNNConfiguration) receivedConf;
		if(conf.descriptors != null) {
			BufferedReader br = getAliasedBufferedReader(filename);
			bw = getAliasedBufferedWriter(filename+".new");
			String line;
			while((line = br.readLine())!=null){
				String tokenes[] = line.split("\\s+");
				bw.write(tokenes[0]);
				for(int j = 0; j < tokenes.length-1;j++) {
					if(j< dtDescriptors.getColumnSize() && !conf.descriptors.containsKey(dtDescriptors.getDescriptorName(j)))continue;
					bw.write("\t"+tokenes[j+1]);
				}
				bw.write("\n");
			}

			bw.close();
			br.close();
			File f = getAliasedFile(filename+".new");
			File dest = getAliasedFile(filename);
			f.renameTo(dest);
		}

		return savedRecords.size();
	}

	@Override
	protected void saveOneRow(int mol, DescriptorsTable dtDescriptors, String[] values, Writer writer) throws IOException {

		writer.append("mol_id" +dtDescriptors.getMoleculeUniqueId(mol)+"_stero"+dtDescriptors.getMoleculeIdStereochem(mol)+VALUESSEPARATOR);

		String descs[]=dtDescriptors.getScaledValuesString(mol);

		for(int i = 0;i<descs.length;i++)
			writer.append(descs[i]+VALUESSEPARATOR);

		for(int i=0; values != null && i <values.length; i++)
			writer.append((values[i] == null ? LabelsTable.MISSED_VALUES: values[i]) + (i !=values.length -1? VALUESSEPARATOR:""));

		writer.append(OSType.endLine());
	}

	public ASNNServer(){
		supportedTaskType = QSPRConstants.ASNN;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}


}
