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
import qspr.metaserver.configurations.DescriptorsMORDREDConfiguration;
import qspr.metaserver.configurations.DescriptorsMORDREDConfiguration.Descr;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

public class MORDREDDescriptorsServer extends DescriptorsAbstractExecutableServer
{

	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration,
			int start,int size) throws Exception
	{
		DescriptorsMORDREDConfiguration c = (DescriptorsMORDREDConfiguration) receivedConfiguration;
		String commands[]={null,"-m", "mordred", datain+".sdf", c.mordred == Descr.D3?"-3":"", "-s", "-o", dataout};
		saveMolecules(dtMolecules, datain+".sdf",QSPRConstants.SDFNOAROM_WITHH,start,size);
		runPython(commands, dataout, CONDA.RDKIT, size <10 ?10:size);
		DataTable res = readStandardOutputResults(getResults(),dataout,false);
		return res;
	}

	@Override
	int getBatchSize(){
		return 10;
	}

	public MORDREDDescriptorsServer()
	{
		supportedTaskType = DescriptorsConfiguration.MORDRED;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		DELIMITER = ",";
		repostSize = 1000;
	}


}
