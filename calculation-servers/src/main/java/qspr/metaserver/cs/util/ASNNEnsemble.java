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

package qspr.metaserver.cs.util;

import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import qspr.metaserver.CalculationServer;
import qspr.metaserver.configurations.ASNNConfiguration;
import qspr.metaserver.configurations.CompressedObject;
import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;

import com.eadmet.utils.NumericalValueStandardizer;

public class ASNNEnsemble extends NeuralEnsemble {

	protected byte[] getRanksByte(int size, float vals[]){
		byte[] ranks=new byte[vals.length];
		int cases=vals.length/size;

		for(int i=0;i<cases;i++){
			Map <Integer,Float> m = new HashMap<Integer,Float>();
			for(int j=0;j<size;j++)
				m.put(j, vals[i*size+j]);	

			Map  <Integer,Float>sorted=sortByValuesDescent(m);

			Set<Integer> set=sorted.keySet();

			Iterator<Integer> it=set.iterator();

			byte j = (byte)(-size/2);

			for(int k=0;k<size;k++,j++){
				int nn=it.next();
				ranks[i*size+nn]=j;
			}
		}

		return ranks;		
	}


	private void addCorrel(double []predictions,byte[] thisRanks, byte[] trainRanks, double biases[],CalculationServer server) {
		int outputs = outputs();
		int samples=thisRanks.length/(ensemble.length*outputs);

		for(int i=0;i<samples;i++){
			if(samples > 10 && i> 0 && (i%(samples/10) == 0))server.setStatus("addCorrel progress: "+(int)(i*100./samples)+"%");
			Map<Integer,Double> correl = getCorrel(thisRanks,trainRanks,i); // getting rank for all outputs simultaneously
			double []cor=getCorrelations(correl, biases);
			for(int j=0;j<outputs;j++)
				predictions[outputs*samples*2+i*outputs+j]=1.-cor[j];  // distances
		}
	}


	private Map<Integer, Double> getCorrel(byte[] thisRanks,
			byte[] trainRanks, int n) {

		Map <Integer,Double> val=new HashMap<Integer,Double>();

		int outputs = outputs();
		int size=ensemble.length*outputs;
		int cases=trainRanks.length/size;

		for(int j=0;j<cases;j++){
			if(j==n && thisRanks==trainRanks)continue;
			double cor=correlation(thisRanks,n*size,trainRanks,j*size,size);
			val.put(j,cor);
		}

		return sortByValuesDescent(val);
	}

	private double correlation(byte[] thisRanks, int startx,
			byte[] trainRanks, int starty,int ensemble) {

		double cor=0,var=0,meanv=(ensemble%2 == 0)?-0.5:0,a,b;


		for(ensemble += startx;startx<ensemble;startx++,starty++){
			a=thisRanks[startx]-meanv;
			b=trainRanks[starty]-meanv;
			cor += a*b;
			var += a*a;
		}

		return cor/var;
	}

	protected void correctASNN(double []predictions,byte[] thisRanks, byte[] trainRanks, double biases[], int k, double parzen, boolean median) {

		int outputs = normY.length/2;
		int samples=thisRanks.length/(ensemble.length*outputs);

		for(int i=0;i<samples;i++){
			Map<Integer,Double> correl = getCorrel(thisRanks,trainRanks,i); // getting rank for all outputs simultaneously

			double []correction= median? getMedian(correl,k,parzen,biases):  getBias(correl,k,parzen,biases);
			for(int j=0;j<outputs;j++)
				predictions[i*outputs+j]+=correction[j];

		}

	}

	public void prepareRanksAndBiases(ASNNConfiguration conf, String fileName) throws IOException {
		LineNumberReader file = new LineNumberReader(new FileReader(fileName));
		readData(file); 
		file = new LineNumberReader(new FileReader(fileName)); // reopening to reset linnumbers

		float newPredictions[]=predict(file); // getting all predictions

		file.close();

		int outputs=normY.length/2;

		byte[] theseRanks=getRanksByte(outputs*ensemble.length, newPredictions);

		for(int i=0;i<newPredictions.length;i++){
			int output=i%outputs;
			newPredictions[i] = (float)((newPredictions[i]-normY[output*2+1])/normY[output*2]);
		}

		double predictions[] = averageAndStd(outputs,ensemble.length,newPredictions);

		conf.modelRanks = new CompressedObject<Object>();
		conf.modelBiases = new CompressedObject<Object>();

		conf.modelRanks.set(theseRanks);
		double biases[] = new double[dataY.length]; // calculate and store errors 

		for(int i=0;i<biases.length;i++)
			biases[i]= dataY[i]>(MISSED_VALUE-1)? MISSED_VALUE:dataY[i]-predictions[i];

			conf.modelBiases.set(biases);
	}

	public void predictASNN(String dataName, String resultFile, ASNNConfiguration conf, CalculationServer server) throws IOException {

		LineNumberReader file = new LineNumberReader(new FileReader(dataName));
		readData(file); 
		file = new LineNumberReader(new FileReader(dataName)); // reopening to reset linnumbers

		server.setStatus("Starting ANN predict ...");

		float vals[]=predict(file); // getting all predictions
		file.close();

		server.setStatus("Finished ANN predict ...");

		BufferedWriter writer = new BufferedWriter(new FileWriter(resultFile));

		int outputs=normY.length/2;

		byte [] trainRanks = null,
				thisRanks = null;

		if(conf.asnn){
			trainRanks = (byte [])conf.modelRanks.get();
			thisRanks =getRanksByte(outputs*ensemble.length, vals);
		}

		for(int i=0;i<vals.length;i++){
			int output=i%outputs;
			vals[i] = (float)((vals[i]-normY[output*2+1])/normY[output*2]);
		}

		server.setStatus("Starting averageAndStd  ...");

		double predictions[] = averageAndStd(outputs,ensemble.length,vals);

		if(conf.asnn && conf.predictionScenario != PredictionScenario.DISTANCE_ONLY){
			server.setStatus("Starting addCorrel  ...");
			addCorrel(predictions,thisRanks,trainRanks,(double[])(conf.modelBiases.get()), server); //  CORREL metric only for the test sets
			server.setStatus("Starting correctASNN  ...");
			correctASNN(predictions, thisRanks, trainRanks,  (double[])(conf.modelBiases.get()), conf.k, conf.parzen, false);  // ASNN correction + calculation of CORREL metric
		}	

		for(int i=0;i<predictions.length/(3*outputs);i++){ // by cases
			writer.write("\nmol_"+i+"\t");

			for(int k=0;k<outputs;k++)
				writer.write("\t0");

			for(int j=0;j<3;j++)
				for(int k=0;k<outputs;k++){
					int nn=j*predictions.length/3 + i*outputs+k;
					writer.write("\t"+NumericalValueStandardizer.getSignificantDigitsStr(predictions[nn],4));
				}

		}
		writer.close();
	}

	public void readData(LineNumberReader file) throws IOException {
		String vals[] = file.readLine().replaceAll("[<>]", "").split("\\s+"); // file name
		int n=Integer.parseInt(vals[0]);
		readDataSet(file, n);
	}

}
