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
import java.util.Arrays;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsMolPrintConfiguration;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;

public class MolPrintServer extends DescriptorsAbstractExecutableServer
{

	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration, int start, int size) throws Exception
	{
		if (receivedConfiguration == null || !(receivedConfiguration instanceof DescriptorsMolPrintConfiguration))
			throw new Exception("Invalid configuration passed, should be instance of MolPrintConfiguration");
		DescriptorsMolPrintConfiguration configuration = (DescriptorsMolPrintConfiguration) receivedConfiguration;

		saveMolecules(dtMolecules, datain, QSPRConstants.MOL2, start, size);

		String commands[]={"python3",getExeFile(),"-i",datain,"-o",dataout,"-d",""+configuration.depth,"-p","mol","-n","1"};
		executeBinary(commands, dataout, 15*size);

		return readResults(dataout,dtMolecules.getRowsSize());
	}

	@Override
	int getBatchSize(){
		return 1000;
	}

	String sort(String feature) {
		String vals[] = feature.split(";");
		Arrays.sort(vals);
		feature = "";
		for(String s: vals)
			feature += feature.length()>0?";"+s:s;
			return feature;
	}

	DataTable readResults(String file, int totalsize) throws IOException {

		BufferedReader in =getAliasedBufferedReader(file);

		String line=null; 
		DataTable res = new DataTable(true);

		while((line=in.readLine())!=null){
			res.addRow();

			String[] features = line.trim().split("\t");
			for (String feature : features) {
				feature = feature.trim();
				feature = sort(feature);
				if(feature.contains("mol"))
					continue;
				res.setValue(feature, 1);
			}
		}
		in.close();

		DataTable dtResults = getResults();
		dtResults.addRowsFrom(res); // using an efficient way to add rows at the end of the table

		dtResults.print(System.out);
		return dtResults;
	}


	public MolPrintServer()
	{
		supportedTaskType = DescriptorsConfiguration.MolPrint;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}

}
