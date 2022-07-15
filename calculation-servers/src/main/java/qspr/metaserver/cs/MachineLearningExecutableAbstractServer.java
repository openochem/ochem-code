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
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.util.List;

import com.eadmet.utils.FileUtils;
import com.eadmet.utils.OCHEMUtils;

import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.util.ExecutableRunner;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

abstract public class MachineLearningExecutableAbstractServer extends MachineLearningAbstractServer{

	public final static String stdout = "stdout";
	protected ExecutableRunner exeRunner;
	protected String DESC = "des";


	/**
	 * To return always card 0 for servers which support card reloading
	 * @param forceCPU
	 * @return
	 */

	protected String getGPUCard(boolean forceCPU) { // always card 0, since environment is used to restrict only one card, which becomes 0
		return "" + (noGPU() || forceCPU ? NO_GPU : gpuCard);
	}

	static boolean noGPU() {
		return OCHEMUtils.findFile("nvidia-smi") == null;
	}

	protected String runPythonWithConda(String [] commands, String resultFile, CONDA env) throws Exception{

		if(getGPUCard(false).equals("0")) {
			commands = OCHEMUtils.insert("CUDA_VISIBLE_DEVICES=" + gpuCard +";",commands);
			commands = OCHEMUtils.insert("export", commands);
			commands = OCHEMUtils.insert("CUDA_DEVICE_ORDER=PCI_BUS_ID;", commands);
			commands = OCHEMUtils.insert("export", commands);
		}
		return exeRunner.runPython(commands, resultFile, env, 0);
	}

	@Override
	void init() throws IOException, InterruptedException{
		exeRunner = new ExecutableRunner(this);
		exeRunner.init();
	}

	@Override
	void cleanup() throws IOException, InterruptedException{
		exeRunner.cleanup();
		exeRunner=null;
	}

	protected String getAliasedPath() {
		return exeRunner.getAliasedPath();
	}

	public void executeBinaryBash(String[] commands, String outputFile) throws IOException, InterruptedException {
		exeRunner.executeBinaryBash(commands, outputFile, null, 0);
	}

	protected File getAliasedFile(String filename) {
		return new File(exeRunner.getAliasedFileName(filename));
	}

	protected String getAliasedFileName(String filename) {
		return exeRunner.getAliasedFileName(filename);
	}

	protected BufferedWriter getAliasedBufferedWriter(String filename) throws IOException {
		return exeRunner.getAliasedBufferedWriter(filename);
	}

	protected BufferedReader getAliasedBufferedReader(String filename) throws IOException {
		return exeRunner.getAliasedBufferedReader(filename);
	}

	protected LineNumberReader getAliasedLineNumberReader(String filename) throws IOException {
		return new LineNumberReader(exeRunner.getAliasedBufferedReader(filename));
	}

	protected void executeBinarySafely(String[] commands, String outputFile) throws IOException, InterruptedException {
		for(int i=0; i< 2; i++)try {
			exeRunner.executeBinary(commands, outputFile);
		}catch(Exception e) {
			Thread.sleep(1000*(i+1));
		}
	}

	protected void executeBinary(String[] commands, String outputFile) throws IOException, InterruptedException {
		exeRunner.executeBinary(commands, outputFile);
	}

	protected void executeBinary(String[] commands, String outputFile, int timeout) throws IOException, InterruptedException {
		exeRunner.executeBinary(commands, outputFile,timeout);
	}

	protected void executeBinary(String[] commands, String outputFile, int timeout,
			PrintWriter out, PrintWriter err) throws IOException, InterruptedException {
		exeRunner.executeBinary(commands, null, outputFile,timeout,out, err);
	}

	protected void executeBinary(String[] commands, String env[], String outputFile, int timeout) throws IOException, InterruptedException {
		exeRunner.executeBinary(commands, env, outputFile, timeout, out, err);
	}

	protected String getExeFile() {
		return  exeRunner.getExeFile();
	}

	protected void setLogStdOut(boolean val) {
		exeRunner.setLogStdOut(val);
	}

	void saveModelToFile(ModelAbstractConfiguration configuration, String aliasedFileName) throws IOException {
		if(configuration.saveModels()) // only if model was saved ...
			FileUtils.saveBytesToFile(configuration.getSavedModelAsBytes(),aliasedFileName);
	}

	DataTable readResultValuesSingle(String filename, List<Integer> optionsNumber, int offset) throws IOException{

		BufferedReader br = exeRunner.getAliasedBufferedReader(filename);
		DataTable dtResult = new DataTable(true);
		dtResult.id = QSPRConstants.PREDICTION_RESULT_COLUMN;

		optionsNumber = optionsNumber == null || optionsNumber.size() == 0? null: optionsNumber;

		String line;
		while ((line = br.readLine()) != null)
			try{
				if(line.contains(QSPRConstants.PREDICTION_RESULT_COLUMN) || line.contains("\"x\"")) continue; // Second is for R---

				dtResult.addRow();
				String[] predictions = line.split(",");

				int classNumber = optionsNumber == null ? 1 :  optionsNumber.get(0);
				Double value = Double.parseDouble(predictions[offset]);
				if(classNumber == 2 ) value = correct(value);  
				dtResult.setValue(QSPRConstants.PREDICTION_RESULT_COLUMN + 0, value);

			}catch(Throwable e) {
				dtResult.getCurrentRow().setError(e.getMessage());
			}

		br.close();
		return dtResult;	
	}


	/**
	 *  Provides rounding of a value to (-0.4; 1.4) interval to match with correct assign of a class
	 * @param x
	 * @return
	 */

	static Double correct(double x) {
		if(x >= 0 && x <= 1) return x;
		if(x > 1) return 0.6 + sigmoid(x -1)*0.8;
		return  0.4 - sigmoid(-x)*0.8;
	}

	static double sigmoid(double x) {
		return 1./( 1 + Math.pow(Math.E,(-1*x))); 
	}


}

