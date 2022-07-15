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

import java.io.BufferedReader;
import java.io.IOException;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;

public class ALogPSServer extends DescriptorsAbstractExecutableServer
{
	final static String dataFile = "data_in.sdf";
	final static String outs = "result.txt";

	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration,int start,int size) throws IOException, InterruptedException 
	{
		saveMolecules(dtMolecules, dataFile,QSPRConstants.SDF,start,size);
		String[] commands = new String[] { getExeFile(), "-f", dataFile, "-o", outs };
		executeBinary(commands, outs,size<30?10:size); // 10 seconds per molecule maximum
		return readResults(size);
	}

	@Override
	protected int getBatchSize(){
		return 1000;
	}

	private DataTable readResults(int moleculesNumber) throws IOException
	{
		DataTable dtResults = getResults();
		boolean first=false;
		String line;

		BufferedReader input = getAliasedBufferedReader(outs);
		input.readLine();
		input.readLine();

		for (int i = 0; i < moleculesNumber; i++)
		{
			line = input.readLine();
			dtResults.addRow();

			if(line==null || line.length()==0)
				throw new IOException("ALOGPS crashed"); // stop any further analysis!

			line=line.toLowerCase();

			first=true;
			if (line.contains("error"))
			{
				dtResults.getCurrentRow().setError(line.substring(6));
				out.println(line);
				continue;
			}
			String[] val = line.split("\\s+");
			if(val.length<3){
				dtResults.getCurrentRow().setError("program did not provide any values");
				continue;
			}
			dtResults.setValue("ALogPS_logP", Double.valueOf(val[1]));
			dtResults.setValue("ALogPS_logS", Double.valueOf(val[2]));
		}

		line=input.readLine(); // empty line
		line=input.readLine(); // Signature line

		input.close();

		if(line == null || line.indexOf("Estimated logP")==-1)
			throw new IOException("Number of processed molecules does not correpond to the submitted one.");

		if(!first){
			dtResults.addRow();
			dtResults.getCurrentRow().setError("program crashed");

		}
		return dtResults;
	}

	public ALogPSServer()
	{
		supportedTaskType = "ALogPS";
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}
}
