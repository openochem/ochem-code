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

import java.io.BufferedWriter;
import java.util.Map;

import qspr.metaserver.configurations.FSMLRConfiguration;
import qspr.metaserver.configurations.LinearConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.utils.QSPRConstants;

public class FSMLRServer extends LinearAbstractServer
{

	@Override
	protected Map<Integer, Double> calculateCoefficients(
			DescriptorsTable dtDescriptors, LabelsTable dtExpValues,
			LinearConfiguration receivedConf) throws Exception {

		BufferedWriter writer = getAliasedBufferedWriter(DATAFILE);
		writeSimpleIgorFormatData(dtDescriptors, dtExpValues, writer);
		writer.close();

		FSMLRConfiguration configuration = (FSMLRConfiguration) receivedConf;

		String[] commands = new String[] { getExeFile(), "-in", DATAFILE, "-model", MODEL, "-mdd",
				configuration.delta + "", "-nfolds", configuration.nfolds + "", "-start", configuration.start + "",
				"-def", configuration.extraFactor + "", "-dnf" + configuration.numFactor,
				"-shr" + configuration.shrinkage };
		executeBinary(commands, OCHEM);

		return getCoefficients(MODEL,true);
	}


	public FSMLRServer()
	{
		supportedTaskType = QSPRConstants.FSMLR;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}

}
