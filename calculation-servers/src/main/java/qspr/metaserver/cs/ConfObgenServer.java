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

import qspr.dao.Various;
import qspr.metaserver.configurations.StructureOptimisationConfiguration;
import qspr.metaserver.util.ExecutableRunner;
import qspr.workflow.utils.QSPRConstants;

public class ConfObgenServer extends ConfBabelServer{

	@Override
	protected String convertedStructure(String molecule, StructureOptimisationConfiguration config) throws Exception {

		String sdf = null;

		try {
			String[] commands = { ExecutableRunner.findExecutable("obgen"), server.getAliasedFileName(sdfFile) };

			FileUtils.saveStringToFile(molecule, server.getAliasedFileName(sdfFile));
			server.executeBinary(commands, ExecutableRunner.stdout, TIMEOUT);
			sdf = FileUtils.getFileAsString(server.getAliasedFileName(ExecutableRunner.stdout));

			if(Various.molecule.compareMolecules(molecule, sdf) == null) return sdf;

			FileUtils.saveStringToFile(Various.molecule.addHydrogensAndRemoveAll(molecule), server.getAliasedFileName(sdfFile));
			server.executeBinary(commands, ExecutableRunner.stdout, TIMEOUT);
			sdf = FileUtils.getFileAsString(server.getAliasedFileName(ExecutableRunner.stdout));
		}catch(Exception e ) {
		}

		if(sdf != null && Various.molecule.compareMolecules(molecule, sdf) == null) return sdf;

		return super.convertedStructure(molecule, config);
	}

	public ConfObgenServer(){
		supportedTaskType = QSPRConstants.OBGEN;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		repostSize=100;
	}

}
