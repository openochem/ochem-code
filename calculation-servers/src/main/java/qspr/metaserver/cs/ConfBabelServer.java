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

import com.eadmet.utils.FileUtils;

import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.metaserver.util.ExecutableRunner;
import qspr.workflow.utils.QSPRConstants;

public class ConfBabelServer extends StructureConversionExecutableAbstractServer
{
	public final static String DATAIN = "i.sdf";
	public final static String DATAOUT = "o.sdf";
	public final static String MEND = "M  END";

	public ConfBabelServer()
	{
		supportedTaskType = QSPRConstants.OBABEL;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		repostSize=500;
	}

	@Override
	String convertedStructure(String molecule, StructureOptimisationConfiguration config) throws Exception {
		String[] commands = {ExecutableRunner.findExecutable("obabel"), "-i" + "sdf", DATAIN, "-o" + "sdf", "-O", server.getAliasedFileName(DATAOUT), "--gen3D", "-h"};
		FileUtils.saveStringToFile(molecule, server.getAliasedFileName(DATAIN));
		String file = null;
		server.runPython(commands, DATAOUT, null, TIMEOUT);
		file = FileUtils.getFileAsString(server.getAliasedFileName(DATAOUT));
		if(!file.contains("nan")) return file;
		throw new IOException("Failed to convert");
	}

}
