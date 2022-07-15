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
import java.io.IOException;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsPaDEL2Configuration;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

public class PaDEL2DescriptorsServer extends DescriptorsAbstractExecutableServer
{
	protected static final String CFG = "config.cfg";
	protected static final String INPUT = "molecules.sdf";

	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration,
			int start,int size) throws Exception
	{
		DescriptorsPaDEL2Configuration conf = (DescriptorsPaDEL2Configuration)receivedConfiguration;
		createCfg(conf);
		saveMolecules(dtMolecules,getAliasedFileName(INPUT),QSPRConstants.SDFNOAROM_WITHH, start, size);
		String[] commands = new String[] {"java","-cp",getAliasedPath()+"*:" + 
				javaClassPath,"padeldescriptor.PaDELDescriptorApp", "-config",CFG};
		exeRunner.executeBinaryBash(commands, dataout, null, size <= 2 ? 120:60*(size+1)); 
		return readStandardOutputResults(getResults(), dataout, false);
	}

	@Override
	int getBatchSize(){
		return 8;
	}

	protected void createCfg(DescriptorsPaDEL2Configuration configuration) throws IOException {

		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		writer.write("Directory="+getAliasedPath());
		writer.write("\nDescriptorFile="+dataout);
		writer.write("\nCompute2D="+configuration.bit(0));
		writer.write("\nCompute3D="+configuration.bit(1));
		writer.write("\nComputeFingerprints="+configuration.bit(2));
		writer.write("\nRemoveSalt=false");
		writer.write("\nDetectAromaticity=true");
		writer.write("\nStandardizeTautomers=false");
		writer.write("\nStandardizeNitro=false");
		writer.write("\nTautomerFile=");
		writer.write("\nRetain3D=true");
		writer.write("\nLog=false");
		writer.write("\nConvert3D=FORCEFIELD_NO");
		writer.write("\nMaxThreads=1");
		writer.write("\nMaxJobsWaiting=200");
		writer.write("\nMaxRunTime=60000"); // im ms
		writer.write("\nMaxCpdPerFile=5000");
		writer.write("\nRetainOrder=true");
		writer.write("\nUseFilenameAsMolName=false");
		writer.write("\n");
		writer.close();
	}


	public PaDEL2DescriptorsServer()
	{
		supportedTaskType = DescriptorsConfiguration.PADEL2;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		repostSize = 32;
		startPosition = 1;
		DELIMITER = ",";
	}

}
