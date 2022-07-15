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

@XmlRootElement(name = "xgboost-configuration")
public class XGBOOSTConfiguration extends ModelAbstractConfiguration{

	private static final long serialVersionUID = 1L;

	public Float eta = 0.2f;

	public Integer depth = 8;

	public Integer rounds = 1024;

	public static final String LINEAR = "reg:linear", LOGISTIC = "reg:logistic", BINARY = "binary:logistic";

	public String objective = LINEAR;

	public Float lambda = 1f;

	public String toString()
	{
		return super.toString() +  " trees: " + rounds + " with max depth of " + depth + " eta=" + eta;
	}

	@Override
	public boolean isSupportRegression() {
		return true;
	}

	@Override
	public String getDefaultName() {
		return QSPRConstants.XGBOOST;
	}
}
