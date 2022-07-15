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

package qspr.metaserver.cs;

import com.eadmet.utils.FileUtils;

import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.workflow.utils.QSPRConstants;

public class ConfBalloonServer extends StructureConversionExecutableAbstractServer
{

	public ConfBalloonServer()
	{
		supportedTaskType = QSPRConstants.BALLOON;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		repostSize=1000;
	}

	@Override
	protected String convertedStructure(String molecule, StructureOptimisationConfiguration config) throws Exception {

		String[] commands = { server.getExeFile(), "-c", "1", "-f", "MMFF94.mff", server.getAliasedFileName(sdfFile), "out.sdf" };

		FileUtils.saveStringToFile(molecule, server.getAliasedFileName(sdfFile));

		server.executeBinaryBash(commands, "out.sdf", null, TIMEOUT);

		return FileUtils.getFileAsString(server.getAliasedFileName("out.sdf"));

	}

}
