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

package qspr.metaserver.configurations;

import java.util.LinkedHashMap;

import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.workflow.utils.QSPRConstants;

@XmlRootElement(name = "asnn-configuration")
public class ASNNConfiguration extends MultiLearningAbstractConfiguration 
{
	private static final long serialVersionUID = 1L;
	public static final int MOMENTUM = 0;
	public static final int SUPERSAB = 1;
	public static final int RPROP = 2;
	public static final int QUICKPROP = 3;
	public static final int DIFF_EQUATION = 4;
	public static final int QUICKPROP_2 = 5;
	public static final int LEVENBERG = 6;

	public Boolean asnn = true;
	public Integer ensemble = 64;
	public Integer neurons = 3;
	public Integer iterations = 1000;
	public Integer training = SUPERSAB;
	public String additionalParam = "PARTITION=3,SELECTION=2";
	public Integer k=0;
	public Double parzen=0.;

	@XmlTransient
	public Integer models=1;

	public Integer libraryOutput; // which output we should use in case if provided model is multi-model	

	public CompressedObject<Object> libraryModel; // the previous "mother" model which will be updated
	public CompressedObject<Object> modelRanks;  // stores ranks
	public CompressedObject<Object> modelBiases; // stores biases
	public LinkedHashMap<String,Integer> descriptors = null;

	public String toString()
	{
		String p=getMethodName()+", "+Math.abs(iterations)+" iterations, "+neurons+" neurons "+" ensemble="+ensemble+ (k>0? " k="+k:"");
		return (additionalParam!=null && additionalParam.length()>0? p+" additional param "+additionalParam:p) +
				super.toString();
	}

	private String getMethodName()
	{
		String res = "";
		switch (training)
		{
		case MOMENTUM: res = "Momentum"; break;
		case SUPERSAB: res = "Supersab"; break;
		case RPROP: res = "Rprop"; break;
		case QUICKPROP: res = "QuickProp"; break;
		case DIFF_EQUATION: res = "Diff. Equation"; break;
		case LEVENBERG: res = "Levenberg"; break;
		}
		return res;
	}

	/**
	 * Creates configuration to perform 1 iteration and thus
	 * calculate correct parameters of the ASNN  
	 */

	public void createASNNApplyModelConfiguration(){
		models=3; // i.e., we will need to read the old and write a new model file
		asnn=true; // asnn model is enabled
		iterations=1; // one iteration to recalculate all arrays and to perform proper ASNN
		k=0;
		parzen=0.;
		additionalParam=null; // no additional params
	}

	@Override
	public boolean isSupportRegression() {
		return true;
	}

	@Override
	public String getDefaultName() {
		return asnn?QSPRConstants.ASNN:QSPRConstants.ANN;
	}

	@Override
	public boolean isSupportPredicates() {
		return true;
	}

	@Override
	public void setIterations(int iterations) {
		this.iterations = iterations;
	}

}

