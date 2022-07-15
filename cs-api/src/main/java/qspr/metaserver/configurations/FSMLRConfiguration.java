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

@XmlRootElement(name = "fsmlra-configuration")
public class FSMLRConfiguration extends LinearConfiguration
{
	private static final long serialVersionUID = 3L;

	public float shrinkage = 1f;
	public float delta = 20f;
	public float start = 0f;
	public float numFactor = 1.0f;
	public float extraFactor = 5.0f;
	public int nfolds = 10;

	@Override
	public String getDefaultName() {
		return QSPRConstants.FSMLR;
	}

}
