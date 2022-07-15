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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;


public class NeuralEnsemble {

	NeuralNetwork ensemble[];
	float normX[],normY[];
	int dataX[]; // contain line number in file 
	float dataY[];
	final public double MISSED_VALUE = 999999.;

	public NeuralEnsemble(){}

	public int outputs(){
		return normY.length/2;
	}

	public void readDataNN(LineNumberReader file,boolean training) throws IOException{
		String vals[] = file.readLine().replaceAll("[<>]", "").split("\\s+"); // file name
		int trainN=0,testN=0;

		for(int i=0;i<vals.length;i++){
			if(i==0)trainN+=Integer.parseInt(vals[i]);
			else
				testN+=Integer.parseInt(vals[i]);;
		}

		readDataSet(file, trainN);
		if(training)return;

		readDataSet(file, testN);

	}


	void readDataSet(LineNumberReader file,int testN) throws IOException{

		dataX = new int[testN];
		dataY = new float[normY.length*testN/2];

		int j=0;
		String line,vals[];

		for(int i=0;i<testN;){

			line=file.readLine().replaceAll("[<>]", "");
			if(line==null)throw new IOException("End of file has reached for test="+testN+" "+i);
			if(line.indexOf("mol_")==-1)continue;
			dataX[i++]=file.getLineNumber();
			vals = line.split("\\s+"); // file name

			if(vals.length == normX.length/2 +1)continue; // only values; no output values

			if(vals.length!=(normX.length+normY.length)/2+1)throw new IOException("The number of data columns was changed found: "+
					(vals.length-1)+" != x: "+normX.length/2+" y: "+normY.length/2+" total="+(normX.length+normY.length)/2);

			int dataLength=normX.length/2+1; // 1 + mol_

			for(int k = 0 ; k < normY.length/2 ;k++)
				dataY[j++]=Float.parseFloat(vals[dataLength+k]);

		}

	}



	protected double[] getCorrelations(Map<Integer, Double> sorted, double biases[]) {

		int outputs=normY.length/2;
		double cor[]=new double[outputs];

		Set<Integer> set=sorted.keySet();

		for(int i=0;i<outputs;i++){
			Iterator<Integer> it=set.iterator();
			while(it.hasNext()){
				Integer n=it.next();
				if(biases[n*outputs+i]==MISSED_VALUE)continue;
				double corr=sorted.get(n);
				if(corr<=0)corr=0;
				cor[i]=corr;
				break;
			}
		}

		return cor;
	}


	protected double[] getBias(Map<Integer,Double> sorted, int k, double parzen, double biases[]) {

		Set<Integer> set=sorted.keySet();

		int outputs=normY.length/2;
		double cor[]=new double[outputs];

		for(int i=0;i<outputs;i++){

			double sumV=0,sumC=0;
			Iterator<Integer> it=set.iterator();

			for(int j=0;j<k && it.hasNext();){
				Integer n=it.next();
				if(biases[n*outputs+i]==MISSED_VALUE)continue;

				double corr=sorted.get(n);
				if(corr<=0)continue;
				corr = Math.exp(-parzen*parzen/(corr*corr));
				sumV += biases[n*outputs+i]*corr;
				sumC += corr;
				j++; 
			}
			cor[i]=sumC>0?sumV/sumC:0;
		}

		return cor;
	}

	protected double[] getMedian(Map<Integer,Double> sorted, int k, double parzen, double biases[]) {

		Set<Integer> set=sorted.keySet();

		int outputs=normY.length/2;

		double cor[]=new double[outputs];
		double val[]=new double[k];

		for(int i=0, j=0;i<outputs;i++){

			Iterator<Integer> it=set.iterator();

			for(j=0;j<k && it.hasNext();){
				Integer n=it.next();
				if(biases[n*outputs+i]==MISSED_VALUE)continue;
				val[j++] = biases[n*outputs+i];
			}
			if(j<3)
				cor[i]= j == 0? 0: j == 1? val[0] : (val[0] + val[1])/2;
				else
				{
					Arrays.sort(val, 0, j);
					cor[i] = j % 2 == 1? val[(j+1)/2] :  (val[j/2] + val[1+j/2])/2;
				}
		}

		return cor;
	}


	protected short[] getRanks(int size, float vals[]){
		short[] ranks=new short[vals.length];
		int cases=vals.length/size;

		for(int i=0;i<cases;i++){
			Map <Integer,Float> m = new HashMap<Integer,Float>();
			for(int j=0;j<size;j++)
				m.put(j, vals[i*size+j]);	

			Map  <Integer,Float>sorted=sortByValuesDescent(m);

			Set<Integer> set=sorted.keySet();

			Iterator<Integer> it=set.iterator();

			for(short j=0;j<size;j++){
				int nn=it.next();
				ranks[i*size+nn]=j;
			}
		}

		return ranks;		
	}


	@SuppressWarnings({ "unchecked", "rawtypes" })
	public static <K extends Comparable,V extends Comparable> Map<K,V> sortByValuesDescent(Map<K,V> map){
		List<Map.Entry<K,V>> entries = new LinkedList<Map.Entry<K,V>>(map.entrySet());

		Collections.sort(entries, new Comparator<Map.Entry<K,V>>() {

			@Override
			public int compare(Entry<K, V> o1, Entry<K, V> o2) {
				return o2.getValue().compareTo(o1.getValue());
			}
		});

		//LinkedHashMap will keep the keys in the order they are inserted
		//which is currently sorted on natural ordering
		Map<K,V> sortedMap = new LinkedHashMap<K,V>();

		for(Map.Entry<K,V> entry: entries){
			sortedMap.put(entry.getKey(), entry.getValue());
		}

		return sortedMap;
	}




	protected double[] averageAndStd(int outputs, int ensemble, float[] vals) {

		int size=vals.length/ensemble;
		double av[]=new double[size*3]; // also reserves place for STD and correl values

		int n=0;

		for(int i=0;i<size/outputs;i++) //samples
			for(int j=0;j<ensemble;j++) //averaging by ensemble
				for(int k=0;k<outputs;k++,n++){ // all outputs are covered
					av[i*outputs+k]+=vals[n];
					av[size+i*outputs+k]+=vals[n]*vals[n];
				}

		for(int j=0;j<size;j++){
			av[j]/=ensemble;
			av[j+size]=Math.sqrt((av[j+size]-av[j]*av[j]*ensemble)/(ensemble-1));

		}

		return av;
	}

	public float[] predict(LineNumberReader file) throws IOException {

		float set[]=new float[normX.length/2];

		float predictions[] = new float[ensemble.length*dataX.length*normY.length/2]; //for all length of predictions
		int r=0;
		for(int j=0;j<dataX.length;j++){
			readOne(dataX[j],set,file);
			for(int i=0;i<ensemble.length;i++){
				float val[]=ensemble[i].predict(set);
				for(int k=0;k<val.length;k++)
					predictions[r++]=(float)val[k];
			}
		}
		return predictions;
	}

	float [] readOne(int linenumber,float set[], LineNumberReader file) throws IOException{
		String line=file.readLine();

		while(file.getLineNumber() != linenumber)
			line=file.readLine();

		line=line.replaceAll("[<>]", "");
		if(line==null)throw new IOException("End of file has reached for line"+linenumber);
		String vals[] = line.split("\\s+"); // file name

		for(int k=0;k<normX.length/2;k++)
			set[k]=Float.parseFloat(vals[k+1])*normX[2*k]+normX[2*k+1];

		return set;
	}

	public void readModelNN(BufferedReader model) throws IOException{
		model.readLine(); // file name
		normX = readNorm(model.readLine()); // normalisation of X-inputs
		normY = readNorm(model.readLine()); // normalisation of Y-values
		String vals[] = model.readLine().split("\\s+"); // file name
		ensemble = new NeuralNetwork[Integer.parseInt(vals[0])];
		int hidden = Integer.parseInt(vals[5]);

		for(int i=0;i<ensemble.length;i++){
			ensemble[i] = new NeuralNetwork(normX.length/2, hidden, normY.length/2);
			ensemble[i].readWij(model);
		}
	}

	float []readNorm(String val)throws IOException{
		String vals[]=val.split("\\s+");
		float array[]= new float[2*Integer.parseInt(vals[0])];
		for(int i=1;i<=array.length;i++)array[i-1]=Float.parseFloat(vals[i]);
		return array;
	}

}
