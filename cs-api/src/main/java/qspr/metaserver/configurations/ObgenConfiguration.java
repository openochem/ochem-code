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

@XmlRootElement(name = "obgen-configuration")
public class ObgenConfiguration extends StructureOptimisationConfiguration
{
	private static final long serialVersionUID = 1L;

	public ObgenConfiguration(boolean bypassCache)
	{
		super(bypassCache);
	}

	public ObgenConfiguration()
	{

	}

	@Override
	public String getTaskType() 
	{
		return QSPRConstants.OBGEN;
	}

	@Override
	public String toString() {
		return "3D by " + QSPRConstants.OBGEN;
	}
}
