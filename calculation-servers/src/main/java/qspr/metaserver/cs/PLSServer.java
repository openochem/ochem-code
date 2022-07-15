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
import java.io.FileReader;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.HashMap;
import java.util.List;

import qspr.metaserver.configurations.LinearConfiguration;
import qspr.metaserver.configurations.PLSConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.util.ExecutableRunner;
import qspr.metaserver.util.ExperimentalDesignHelper;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.FileUtils;

public class PLSServer extends LinearAbstractServer {

	final static String OptimiseLatentComponents = "R2.R";
	final static String GetWeights = "getWeights.R";
	final static String OUTPUT = "R2.csv";
	final static String COEFICIENTS = "coef.csv";
	final static String WEIGHTS = "weights.csv";

	public PLSServer() 
	{
		supportedTaskType = QSPRConstants.PLS;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}


	@Override
	public boolean isCritical(String message) {
		if(message.contains("More segments than observations requested"))
			throw new CriticalException("R failed to process these data. Try to use logarithm of the property (select log unit when starting a task) or decrease the number of used descriptors and/or to use larger dataset.");
		return false;
	}

	@Override
	protected HashMap<Integer, Double> calculateCoefficients(
			DescriptorsTable dtDescriptors, LabelsTable dtExpValues,
			LinearConfiguration configuration) throws Exception {


		PLSConfiguration plsConf = (PLSConfiguration) configuration;

		String method = "simpls";

		out.println("A: " + plsConf.numLatentVariables);

		saveTrainingSetDataforR(DATAFILE,dtDescriptors, dtExpValues); 

		if (plsConf.numLatentVariables == 0) {
			FileUtils.saveStringToFile(ExperimentalDesignHelper.scriptOptimiseNumberOfLatentVariables(getAliasedFileName(DATAFILE), getAliasedFileName(OUTPUT), method),
					getAliasedFileName(OptimiseLatentComponents));

			String commands[]={ExecutableRunner.findExecutable("R"),"--slave","-f",OptimiseLatentComponents};
			runPythonWithConda(commands, OUTPUT, null);

			out.println(getAliasedFileName(OUTPUT));
			plsConf.numLatentVariables = getNumberOfOptimalVariables(getAliasedFileName(OUTPUT),25); // maximum of PLS 25 components will be allowed	
		}

		out.println("B: " + plsConf.numLatentVariables);

		FileUtils.saveStringToFile(ExperimentalDesignHelper.scriptWeights(getAliasedFileName(DATAFILE), getAliasedFileName(WEIGHTS), getAliasedFileName(COEFICIENTS), plsConf.numLatentVariables, method),
				getAliasedFileName(GetWeights));

		String commands[]={ExecutableRunner.findExecutable("R"),"--slave","-f",GetWeights};

		for(int step=0 ; step < 1000; step++){
			try{
				runPythonWithConda(commands, WEIGHTS, null);
				break;
			}catch(Throwable e){
				String message="cannot allocate vector of size";
				// too large dataset -- decreasing size to get one fitting in the memory

				setStatus("catched: "+e.getLocalizedMessage());

				int index = e.getMessage().indexOf(message);

				if(index != -1){
					String vals[]=e.getMessage().substring(index+message.length()+1).trim().split("\\s+");
					if(vals.length>1)setStatus("values: "+vals[0]+" : "+vals[1]);
					if(vals.length>1 && vals[1].equalsIgnoreCase("GB") && Double.parseDouble(vals[0])>1.){
						RandomAccessFile f = new RandomAccessFile(getAliasedFileName(DATAFILE), "rw");
						setStatus("Decreasing file by 5% to fit in the memory, step=" + step + " size="+f.length()+" for error "+message+" "+vals[0]+" "+vals[1]);
						long length = (long)(f.length() *0.95);  // decreasing size by 5% to avoid the error
						do {                     
							length -= 1;
							f.seek(length);
						} while(f.readByte() != 10 && length > 0);
						f.setLength(length+1);
						f.close();
						continue;
					}
				}

				throw new UserFriendlyException(e.getMessage());
			}
		}

		//loop
		double[][] tempMatrix = ExperimentalDesignHelper.readCSVFile(getAliasedFileName(COEFICIENTS), 1, plsConf.numLatentVariables, false, ",");

		for (int i = 0; i < tempMatrix.length; i++)
			for (int j = tempMatrix[i].length - 1; j > 0; j--)
				tempMatrix[i][j] -= tempMatrix[i][j-1];

		plsConf.lvMatrix = tempMatrix;

		List<Double> coefficients = ExperimentalDesignHelper.coefficientList(getAliasedFileName(WEIGHTS));
		HashMap<Integer, Double> coeff =  new HashMap<Integer, Double>();

		for(int i = 0;i<coefficients.size();i++) // we store only non-zero coefficients
			if(i==0 || coefficients.get(i)!=0)coeff.put(i, coefficients.get(i));

		return coeff;
	}

	private int getNumberOfOptimalVariables(String fileLocation, int maximumPLSComponents) throws IOException{

		BufferedReader br = null;
		Integer plsComponents = null;
		Double rmse = Double.POSITIVE_INFINITY;

		try {
			br = new BufferedReader(new FileReader(fileLocation));

			br.readLine();
			String line = br.readLine();
			String temp[] = line.split("\\s+");

			for (int i = 2; i < maximumPLSComponents && i< temp.length; i++){
				double val = Double.parseDouble(temp[i]);
				if ( rmse > 1.04*val ) { //??? what is the role of magic number 1.04 ???
					rmse = val;
					plsComponents = i-1;
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		}finally{
			if(br!=null)br.close();
		}

		if(plsComponents==null)throw new UserFriendlyException("Could not read the number of PLS components.");

		return plsComponents;
	}


}
