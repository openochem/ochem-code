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

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "skl-configuration")
public class SKLConfiguration extends ModelAbstractConfiguration{
	private static final long serialVersionUID = 1L;

	public SKLMethod method;

	public enum SKLMethod
	{
		ALL_CV, ALL, ADA_BOOST, ADA_BOOST_TREES, ARD, BAGGING, BAGGING_TREES, BAYESIAN_RIDGE, ELASTIC_NET_CV, ELASTIC_NET, EXTRA_TREES, 
		GAUSSIAN_PROCESS, GRADIENT_BOOSTING, HUBER, KERNEL_RIDGE, K_NEIGHBORS, LASSO_LARS_CV, LASSO_LARS_IC, LINEAR, LOGISTIC_CV, 
		ORTHOGONAL_MATCHING_PURSUIT_CV, RANDOM_FOREST, RANSAC_TREES, RIDGE_CV, SVM, CAT_BOOST, 
		KERNEL_RIDGE_HYPER, GBM
		
		;

		public static boolean isRegressionMethod(){
			return true;
		}

		public static boolean isClassificationMethod(SKLMethod method){
			return  false;
		}
	}

	@Override
	public boolean isSupportRegression(){
		return true;
	}

	public boolean isClassification(){
		return method == SKLMethod.CAT_BOOST;
	}

	@Override
	public String toString(){
		return getDefaultName() + super.toString();
	}

	@Override
	public String getDefaultName() {
		return "SCI:"+method.toString();
	}

}
