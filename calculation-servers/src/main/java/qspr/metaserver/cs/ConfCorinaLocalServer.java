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

import java.io.IOException;

import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.FileUtils;

public class ConfCorinaLocalServer extends StructureConversionExecutableAbstractServer 
{
	@Override
	String convertedStructure(String molecule,
			StructureOptimisationConfiguration config) throws Exception {

		String datafile = "data.sdf", logFile = "log.txt", output = "output.sdf";

		FileUtils.saveStringToFile(molecule, server.getAliasedFileName(datafile));
		String commands[] = { server.getExeFile(), "-i", "t=sdf", "-t", "tracefile=" + logFile, "-d", "wh", datafile, output };

		server.executeBinary(commands, output);
		String mol = FileUtils.getFileAsString(server.getAliasedFileName(output));
		if (mol.length() < 5)
			throw new IOException(supportedTaskType + " failed to process this strcuture");

		return mol;
	}

	public ConfCorinaLocalServer()
	{
		supportedTaskType = QSPRConstants.CORINA;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		repostSize = 5000;
	}

}
