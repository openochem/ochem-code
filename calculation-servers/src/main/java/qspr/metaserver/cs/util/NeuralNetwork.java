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

public class NeuralNetwork {
	/*
	 * Object to store a single Neural Network with one hidden layer, bias and sigmoid activation function
	 */
	
	final static double LIMIT=50.;
	
	public double inHid[], hidOut[];
	float hidden[], outputs[]; // weights of the network, including bias
	
	
	public NeuralNetwork(int inputs, int hidden, int outputs) {
		inHid  = new double[(inputs+1)*hidden];
		this.hidden  = new float[hidden];
		hidOut = new double[(hidden+1)*outputs];
		this.outputs = new float[outputs];
	}
	
	public static void readVal(double []array,String val)throws IOException{
		String vals[]=val.split("\\s+");
		if(Integer.parseInt(vals[0])!=array.length)throw new IOException("Sizes are different ones: "+vals[0]+"\t"+array.length);
		for(int i=1;i<=array.length;i++)array[i-1]=Double.parseDouble(vals[i]);
	}
	
	public void readWij(BufferedReader data) throws IOException{
		String values=data.readLine();
		readVal(inHid,values);
		values=data.readLine();
		readVal(hidOut,values);
	}
	
	float [] predict(float inputs[]){
		forward(inputs,inHid,hidden);
		float val[] = forward(hidden,hidOut,outputs);
		return val;
	}
	
	float [] forward(float val[],double weights[],float res[]){
		int w=0,wij=weights.length/res.length; // number of weights in the layer
		for(int i=0;i<res.length;i++){
			res[i]=0.f;
			for(int j=0;j<wij-1;j++){
				res[i]+=val[j]*weights[w++];
			}
			res[i]+=weights[w++]; // bias value			
			res[i]=(float)(res[i]>LIMIT?1:res[i]<-LIMIT?0:1./(1.+Math.exp(-(res[i]))));
		}
	
		return res;
	}
	
	
}
