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

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import javax.xml.bind.annotation.XmlRootElement;

import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.NumericalValueStandardizer;

@XmlRootElement(name = "libsvm-configuration")
public class LibSvmConfiguration extends ModelAbstractConfiguration
{
	private static final long serialVersionUID = 1L;
	public Boolean useNu = false;
	public Boolean oneClass = false;
	public String oneClassLabel = null;

	public static final int LinearKernel = 0;
	public static final int PolynomialKernel = 1;
	public static final int RBFKernel = 2;
	public static final int SigmoidKernel = 3;

	public Integer svm_type = -1; //SVM Type (0 = C-SVC, 1 = nu-SVC, 3 = epsilon-SVR, 4 = nu-SVR), classification-regression autodetectable
	public Integer kernel_type = RBFKernel;
	public Integer degree = 3;
	public Double gamma = null; // default: 1/num_features
	public Integer coef0 = 0;
	public Double nu = 0.001;
	public Double cost = 1.;
	public Double epsilon = 0.001;
	public Double svrEpsilon = 0.001;
	public Double classWeightRatio = 1.;

	public Double costMax = 10., costMin = -10., costStep = 2.;
	public Double gammaMax = 10., gammaMin = -10., gammaStep = 2.;
	public Double classWeightRatioMax = 1.0, classWeightRatioMin = 0.1, classWeightRatioStep = 0.1;
	public Double svrEpsilonMax = 10., svrEpsilonMin = -16., svrEpsilonStep = 2.;

	public Boolean useWeighting = false;
	public Boolean gridSearch = true;

	public Double gridSearchSetSize = 0.1; //Fraction of training set used in grid search CV

	public Integer PARALLEL; // If parallel calculations should be performed
	public Integer stepStart, stepStop; // start and stop steps
	public Double  bestAccuracy;

	private Double costCurrent = 0.;
	private Double gammaCurrent = 0.;
	private Double svrEpsilonCurrent = 0.;
	private Double classWeightRatioCurrent = 0.;

	private Double costBest = 0.;
	private Double gammaBest = 0.;
	private Double svrEpsilonBest = 0.;
	private Double classWeightRatioBest = 0.;

	private List<Double> classWeightRatioSteps = new ArrayList<Double>();

	public LibSvmConfiguration() {
		scaleTypeX=ScalingType.RANGE;
	}

	public void setType(Integer val){
	}

	public Integer getType(){
		return (useNu == null || !useNu?0:1) + (oneClass == null || !oneClass?0:1);
	}

	public String getActualParametersString()
	{
		String cost = NumericalValueStandardizer.formattedFloatValuesForWEKA(Math.pow(2, costCurrent));
		String gamma = NumericalValueStandardizer.formattedFloatValuesForWEKA(Math.pow(2, gammaCurrent));
		String svrEpsilon = NumericalValueStandardizer.formattedFloatValuesForWEKA(Math.pow(2, svrEpsilonCurrent));
		if (areClassificationData())
			return "c = " + cost + ", gamma = " + gamma + ((useWeighting != null && useWeighting) ? ", class weight ratio = " + NumericalValueStandardizer.getSignificantDigitsStr(classWeightRatio,3) : "");
		else
			return "c = " + cost + ", gamma = " + gamma + ", epsilon = " + svrEpsilon;
	}

	public void copy(LibSvmConfiguration c){
		useNu = c.useNu;
		oneClass = c.oneClass;
		oneClassLabel = c.oneClassLabel;
		svm_type = c.svm_type;
		kernel_type = c.kernel_type;
		degree = c.degree;
		gamma = c.gamma;
		coef0 = c.coef0;
		nu = c.nu;
		cost = c.cost;
		epsilon = c.epsilon;
		svrEpsilon = c.svrEpsilon;
		seed = c.getSeed();
		classWeightRatio = c.classWeightRatio;

		costMax = c.costMax; costMin = c.costMin; costStep = c.costStep;
		gammaMax = c.gammaMax; gammaMin = c.gammaMin; gammaStep = c.gammaStep;
		classWeightRatioMax = c.classWeightRatioMax; classWeightRatioMin = c.classWeightRatioMin; classWeightRatioStep = c.classWeightRatioStep;
		svrEpsilonMax = c.svrEpsilonMax; svrEpsilonMin = c.svrEpsilonMin; svrEpsilonStep = c.svrEpsilonStep;

		useWeighting = c.useWeighting;
		gridSearch = c.gridSearch;

		gridSearchSetSize = c.gridSearchSetSize;
		PARALLEL = c.PARALLEL;
		stepStart = c.stepStart; stepStop = c.stepStop;
		bestAccuracy = c.bestAccuracy;

		costCurrent = c.costCurrent;
		gammaCurrent = c.gammaCurrent;
		svrEpsilonCurrent = c.svrEpsilonCurrent;
		classWeightRatioCurrent = c.classWeightRatioCurrent;

		costBest = c.costBest;
		gammaBest = c.gammaBest;
		svrEpsilonBest = c.svrEpsilonBest;
		classWeightRatioBest = c.classWeightRatioBest;

		classWeightRatioSteps = c.classWeightRatioSteps;

	}

	private void getActualParameters()
	{
		cost = Math.pow(2, costCurrent);
		gamma = Math.pow(2, gammaCurrent);
		svrEpsilon = Math.pow(2, svrEpsilonCurrent);
		classWeightRatio = classWeightRatioCurrent;
	}

	public void saveBestParameters()
	{
		costBest = costCurrent;
		gammaBest = gammaCurrent;
		svrEpsilonBest = svrEpsilonCurrent;
		classWeightRatioBest = classWeightRatioCurrent;
	}

	public void restoreBestParameters()
	{
		costCurrent = costBest;
		gammaCurrent = gammaBest;
		svrEpsilonCurrent = svrEpsilonBest;
		classWeightRatioCurrent = classWeightRatioBest;
		getActualParameters();
	}

	private void resetGridSearch()
	{
		costCurrent = costMin;
		gammaCurrent = gammaMin;
		svrEpsilonCurrent = svrEpsilonMin;
		classWeightRatioCurrent = classWeightRatioMin;

		resetClassWeightSearch();
	}

	private void resetClassWeightSearch()
	{
		// Fancy double-interval search, e.g. [0.1, 0.2], [5, 10]
		classWeightRatioSteps.clear();

		if (useWeighting == null || !useWeighting)
		{
			classWeightRatioSteps.add(1D);
			return;
		}

		double tmp = classWeightRatioMin;
		while (tmp <= classWeightRatioMax)
		{
			classWeightRatioSteps.add(tmp);
			tmp += classWeightRatioStep;
		}

		tmp = 1 / classWeightRatioMax;
		while (tmp <= 1 / classWeightRatioMin)
		{
			if (Math.abs(tmp - 1.0) > 0.001)
				classWeightRatioSteps.add(tmp);
			tmp += classWeightRatioStep / (classWeightRatioMin * classWeightRatioMax);
		}

		Collections.sort(classWeightRatioSteps);
		classWeightRatioSteps.remove(0);
	}

	public int totalGridSteps()
	{
		int steps = 0;
		resetGridSearch();
		while (nextGridStep())
			steps++;
		return steps;
	}

	public void configurationForGridStep(int step){
		int currentStep=0;
		resetGridSearch();
		while (nextGridStep() && currentStep<step)
			currentStep++;
	}

	private boolean nextGridStepRegression()
	{
		costCurrent += costStep;

		if (costCurrent > costMax)
		{
			costCurrent = costMin;
			gammaCurrent += gammaStep;
		}

		if (gammaCurrent > gammaMax)
		{
			gammaCurrent = gammaMin;
			svrEpsilonCurrent += svrEpsilonStep;
		}

		getActualParameters();
		return svrEpsilonCurrent <= svrEpsilonMax;
	}

	private boolean nextGridStepClassification()
	{
		costCurrent += costStep;

		if (costCurrent > costMax)
		{
			costCurrent = costMin;
			gammaCurrent += gammaStep;
		}

		if (gammaCurrent > gammaMax)
		{
			gammaCurrent = gammaMin;
			classWeightRatioCurrent = classWeightRatioSteps.remove(0);
		}

		getActualParameters();
		return (!classWeightRatioSteps.isEmpty());
	}

	private boolean nextGridStep()
	{
		if (svm_type == 0 || svm_type == 1 || svm_type == 2)
			return nextGridStepClassification();
		else
			return nextGridStepRegression();
	}

	public List<Integer> prepareSteps(){
		ArrayList<Integer> ar= new ArrayList<Integer>();
		int n=totalGridSteps();
		for(int i=0;i<n;i++)
			ar.add(i);
		Collections.shuffle(ar,new Random(getSeed()));
		return ar;
	} 


	public String toString()
	{
		if(areMultiLearningData())
			return "Multiple configurations were used." + super.toString();

		if (!isTrainingConfiguration())
		{
			String str = getSVMType() + ", " + getKernelType();

			String gammaStr = "0";
			if (gamma != null)
				gammaStr = NumericalValueStandardizer.formattedFloatValuesForWEKA(gamma);

			if (areClassificationData())
			{
				if (useWeighting != null && useWeighting)
					str += String.format("(C=%s, gamma=%s, class weight=%s)", 
							NumericalValueStandardizer.formattedFloatValuesForWEKA(cost),
							gammaStr,
							NumericalValueStandardizer.formattedFloatValuesForWEKA(classWeightRatio));
				else
					str += String.format("(C=%s, gamma=%s)", 
							NumericalValueStandardizer.formattedFloatValuesForWEKA(cost),
							gammaStr);					
			}		
			else
				str += String.format("(C=%s, gamma=%s, epsilon=%s)", NumericalValueStandardizer.formattedFloatValuesForWEKA(cost),
						NumericalValueStandardizer.formattedFloatValuesForWEKA(gamma), NumericalValueStandardizer.formattedFloatValuesForWEKA(svrEpsilon));
			return str + super.toString();
		}
		else
			return (oneClass ? "one class" : (useNu ? "new algorithm" : "classic algorithm")) + ", " + getKernelType() + super.toString();
	}

	public String getKernelType()
	{
		String kType = "";
		switch (kernel_type)
		{
		case 0:
			kType = "linear";
			break;
		case 1:
			kType = "polynomial";
			break;
		case 2:
			kType = "RBF";
			break;
		case 3:
			kType = "sigmoid";
			break;
		}
		return kType;
	}

	public String getSVMType()
	{
		String svmType = "";
		switch (svm_type)
		{
		case 0:
			svmType = "C-SVC";
			break;
		case 1:
			svmType = "nu-SVC";
			break;
		case 2:
			svmType = "One-class SVM";			
			break;
		case 3:
			svmType = "epsilon-SVR";
			break;
		case 4:
			svmType = "nu-SVR";
			break;
		}
		return svmType;
	}

	@Override
	public boolean isSupportRegression() {
		return true;
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.LIBSVM;
	}

}
