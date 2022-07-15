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

package qspr.metaserver.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

public class ExperimentalDesignHelper {

	private static transient final Logger logger = LogManager.getLogger(ExperimentalDesignHelper.class);

	public static String getRScript(String plsTrainingsFile, String plsApplicationFile, String latentVariablesFile) {

		String rContent = "";

		// I added initialization of seed number to have reproducibility of PLS models too

		rContent += "library(pls)\n";
		rContent += "set.seed("+QSPRConstants.SEED+")\n";
		rContent += "dataFile <- read.table(file=\"" + plsTrainingsFile + "\", header=TRUE, sep=\"\\t\")\n";
		rContent += "dataContent = list(a = as.matrix(dataFile[,dim(dataFile)[2] ]), b = as.matrix(dataFile[1:dim(dataFile)[2]-1]))\n";
		rContent += "plsModel = plsr(dataContent$a ~ dataContent$b, data = dataContent)\n";
		rContent += "\n";
		rContent += "dataFile <- read.table(file=\"" + plsApplicationFile + "\", header=TRUE, sep=\"\\t\")\n";
		rContent += "dataContent = list(a = as.matrix(dataFile[,dim(dataFile)[2] ]), b = as.matrix(dataFile[1:dim(dataFile)[2]-1]))\n";
		rContent += "\n";
		rContent += "outputData = predict(plsModel,newdata = dataContent,type = \"scores\")\n";
		rContent += "write.table(outputData, file = \"" + latentVariablesFile + "\")\n";

		return rContent;
	}


	public static String scriptWeights(String plsTrainingsFile, String resultFile, String coefFile, int latentVariablesNr, String method) {

		String rContent = "";

		// I added initialization of seed number to have reproducibility of PLS models too

		rContent += "library(pls)\n";
		rContent += "set.seed("+QSPRConstants.SEED+")\n";
		rContent += "dataFile <- read.table(file=\"" + plsTrainingsFile + "\", header=TRUE, sep=\"\\t\")\n";
		rContent += "dataContent = list(a = as.matrix(dataFile[,dim(dataFile)[2] ]), b = as.matrix(dataFile[1:dim(dataFile)[2]-1]))\n";
		rContent += "plsModel = mvr(ncomp= " + latentVariablesNr + ", method=\""+method+"\",  dataContent$a ~ dataContent$b, data = dataContent)\n";
		rContent += "\n";
		rContent += "write.table(coef(plsModel, intercept = TRUE), file = \"" + resultFile + "\")\n";
		rContent += "write.table(plsModel$coefficients,file = \"" + coefFile + "\",sep=\",\")\n";

		return rContent;
	}

	public static String scriptOptimiseNumberOfLatentVariables(String plsTrainingsFile, String resultFile, String method) { // to select the number of components

		String rContent = "";

		// I added initialization of seed number to have reproducibility of PLS models too; probably it is required only at this place

		rContent += "library(pls)\n";
		rContent += "set.seed("+QSPRConstants.SEED+")\n";
		rContent += "dataFile <- read.table(file=\"" + plsTrainingsFile + "\", header=TRUE, sep=\"\\t\")\n";
		rContent += "dataContent = list(a = as.matrix(dataFile[,dim(dataFile)[2] ]), b = as.matrix(dataFile[1:dim(dataFile)[2]-1]))\n";
		rContent += "dataFile <- read.table(file=\"" + plsTrainingsFile + "\", header=TRUE, sep=\"\\t\")\n";
		rContent += "rm(dataFile)\n"; // to release memory for the second copy of the dataset, which is not required anymore
		rContent += "plsModel = mvr(dataContent$a ~ dataContent$b, data = dataContent, method=\""+method+"\" , validation=\"CV\", segments=5)\n";
		rContent += "\n";
		rContent += "write.table(RMSEP(plsModel)$val, file = \"" + resultFile + "\")";

		return rContent;
	}

	public static HashMap<Integer, String> dataTableToMatrix(DataTable dt, boolean property) throws IOException {

		HashMap<Integer, String> lines = new HashMap<Integer, String>();

		if (property) {
			lines.put(-1, "");
		}
		else {
			String temp = "";
			for (int i = 0; i < dt.getColumnsSize(); i++) {
				temp += "Desc" + i + "\t";
			}
			lines.put(-1, temp.substring(0, temp.length() -1));
		}
		int count = 0;

		dt.getColumnsSize();

		dt.reset();
		while (dt.nextRow()) {
			String temp = "";
			for (int i = 0; i < dt.getColumnsSize(); i++) {
				temp += dt.getValue(i) + "\t";

			}
			lines.put(count, temp.substring(0, temp.length() -1));
			count++;
		}

		return lines;
	}

	public static double[][] dataTableToArray(DataTable dt) {

		double array[][] = new double[dt.getRowsSize()][dt.getColumnsSize()];

		int count = 0;
		dt.reset();
		while (dt.nextRow()) {
			for (int i = 0; i < dt.getColumnsSize(); i++) {
				array[count][i] = (Double) dt.getValue(i);
			}
			count++;
		}

		return array;
	}

	public static void printList(HashMap<Integer, String> list, String fileLocation, boolean header) {

		BufferedWriter bw;
		try {
			bw = new BufferedWriter(new FileWriter(fileLocation));

			if (header) {
				bw.write(list.get(-1) + "\n");
			}
			for (int i = 0; i < list.size()- 1; i++) {
				bw.write(list.get(i) + "\n");
			}
			bw.flush();
			bw.close();

		} catch (IOException e) {
			e.printStackTrace();
		}
	}


	public static void pca(String infile, String outfile, int principalPropertiesNumber) {
		throw new UserFriendlyException("PCA is not available anymore as part of Experimental Design Server. It should be implemented as new method.");
	}

	/*		
		try {
			DataSource source = new DataSource(infile);

			Instances data = source.getDataSet();

			PrincipalComponents pca = new PrincipalComponents();
			pca.setMaximumAttributeNames(5);
			pca.setMaximumAttributes(principalPropertiesNumber);
			pca.setInputFormat(data);
			Instances newData = Filter.useFilter(data, pca);

			CSVSaver csv = new CSVSaver();
			csv.setInstances(newData);
			csv.setFile(new File(outfile));
			csv.writeBatch();
		}
		catch (Exception e) {
			e.printStackTrace();
		}

	}
	 */	
	public static List<Double> coefficientList(String fileLocation) {

		BufferedReader br;
		List<Double> coefficients = new ArrayList<Double>();

		try {
			br = new BufferedReader(new FileReader(fileLocation));
			String line = null;

			while ((line = br.readLine()) != null) {

				String temp[] = line.split("\" ");

				if (temp.length == 2) {
					coefficients.add(Double.parseDouble(temp[1]));
				}
			}
			br.close();
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		return coefficients;

	}

	public static double[][] readCSVFile(String fileLocation, int initColumn, int numberOfColumns, boolean readFirstLine, String separator) {

		BufferedReader br;
		HashMap<Integer, String[]> line2Vals = new HashMap<Integer, String[]>();
		int init = initColumn;
		int end = initColumn + numberOfColumns;
		int max = 0;

		try {
			br = new BufferedReader(new FileReader(fileLocation));
			String line = null;

			int lineNumber = -1;
			if (readFirstLine) {
				lineNumber++;
			}

			while ((line = br.readLine()) != null) {

				if (lineNumber > -1) {
					String temp[] = line.split(separator);
					max = temp.length;
					line2Vals.put(lineNumber, temp);
				}

				lineNumber++;
			}
			br.close();
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}

		if (max < end || numberOfColumns == 0) {
			end = max;
		}
		if (init > end) {
			init = 0;
			end = 0;
		}

		double matrix[][] = new double[line2Vals.size()][end-init];

		for (int id : line2Vals.keySet()) {
			for (int i = init; i < end; i++) {
				matrix[id][i - init] = Double.parseDouble(line2Vals.get(id)[i]);
			}
		}

		return matrix;

	}

	public static double[][]  calculateDistanceMatrix(double[][] descriptorMatrix) {

		double distanceMatrix[][] = new double[descriptorMatrix.length][descriptorMatrix.length];

		for (int i = 0; i < descriptorMatrix.length; i++) {
			for (int j = i; j < descriptorMatrix.length; j++) {
				double distance = 0;

				for (int n = 0; n < descriptorMatrix[0].length; n++) {
					distance = distance + (descriptorMatrix[i][n] - descriptorMatrix[j][n]) * (descriptorMatrix[i][n] - descriptorMatrix[j][n]);
				}

				distance = Math.sqrt(distance);

				distanceMatrix[i][j] = distance;
				distanceMatrix[j][i] = distance;
			}
		}

		return distanceMatrix;
	}

	public static double findOptimalExponent(double[][] distanceMatrix, int compoundsToselect) 
	{
		logger.info("IN");

		Integer center = 0;
		double score = 10000000;

		for (Integer i = 0; i < distanceMatrix.length; i++) {
			double tempScore = 0;
			for (Integer j = 0; j < distanceMatrix.length; j++) {
				tempScore += distanceMatrix[i][j];
			}

			if (tempScore < score) {
				center = i;
				score = tempScore;
			}
		}

		double exponent = 1;

		int countSolved = distanceMatrix.length;

		while (countSolved > distanceMatrix.length / compoundsToselect && exponent < 500) {
			logger.info(exponent);
			countSolved = 0;
			for (Integer i = 0; i < distanceMatrix.length; i++) {
				if (Math.pow(1 - distanceMatrix[i][center], exponent) > 0.75) {
					countSolved++;
				}
			}
			logger.info("\tc:\t" + countSolved);
			exponent++;
		}

		while (countSolved < distanceMatrix.length / compoundsToselect) {
			logger.info(exponent);
			countSolved = 0;
			for (Integer i = 0; i < distanceMatrix.length; i++) {
				if (Math.pow(1 - distanceMatrix[i][center], exponent) > 0.75) {
					countSolved++;
				}
			}
			logger.info("\tc:\t" + countSolved);
			exponent = exponent -0.1;
		}

		logger.info(countSolved);
		logger.info("OUT");

		return exponent;

	}

	public static double[] getColumn(double[][] matrix, int column) {

		double col[] = new double[matrix.length];

		for (int i = 0; i < matrix.length; i++) {
			col[i] = matrix[i][column];
		}

		return col;

	}

	public static boolean isEven(int number) {
		boolean even = false;
		int temp = number / 2;
		temp = temp * 2;
		if (temp == number) {
			even = true;

		}
		return even;
	}

	public static double normalize(double value, double minVal, double maxVal) {

		double result = (value - minVal) / (maxVal - minVal) ;
		return result;
	}

}
