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

import qspr.workflow.utils.QSPRConstants;

@XmlRootElement(name = "knn-configuration")
public class KNNConfiguration  extends ModelAbstractDataDrivenConfiguration
{
	private static final long serialVersionUID = 2L;
	public static final int EUCLIDIAN = 0;
	public static final int PEARSON = 1;

	public Integer distance = EUCLIDIAN;
	public Integer maxKnn = 100;
	public Integer knn=0; 

	public String version;
	public String cfg;

	public String toString()
	{
		String res = "";
		if (knn != null && knn > 0)
			res += "Fixed KNN = "+knn;
		else
			res += "Max KNN = "+maxKnn;

		if (distance.equals(EUCLIDIAN))
			res += ", Euclidian distance";
		else
			res += ", Pearson distance";

		return res + super.toString();
	}

	@Override
	public boolean isSupportRegression() {
		return true;
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.KNN;
	}
}
