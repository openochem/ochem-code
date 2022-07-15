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

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsMOLD2Configuration;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;

public class MOLD2DescriptorsServer extends DescriptorsAbstractExecutableServer
{

	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration, int start, int size) throws Exception
	{
		if (receivedConfiguration == null || !(receivedConfiguration instanceof DescriptorsMOLD2Configuration))
			throw new Exception("Invalid configuration passed, should be instance of DescriptorsMold2Configuration");
		saveMolecules(dtMolecules, datain, QSPRConstants.SDFNOH, start, size);
		String commands[]={getExeFile(),"-i",getAliasedFileName(datain),"-o",dataout,"<",datain};
		executeBinaryBash(commands, dataout, 15*size);
		return readStandardOutputResults(getResults(),dataout,false);
	}

	@Override
	int getBatchSize(){
		return 100;
	}

	public MOLD2DescriptorsServer()
	{
		supportedTaskType = DescriptorsConfiguration.MOLD2;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		DELIMITER = "\\s+";
		startPosition = 1;
		repostSize = 1000;
	}

}
