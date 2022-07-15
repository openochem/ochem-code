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
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.util.ExecutableRunner;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;

public class PyDescriptorsServer extends DescriptorsAbstractExecutableServer{

	final static String MOLFILE = "molecule.mol2"; 
	final static String OUTPUFILE = "PyDescriptor.csv"; 

	@Override
	protected
	DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration configuration, int start, int batchSize)
			throws Exception {

		saveMolecules(dtMolecules, MOLFILE, QSPRConstants.MOL2, start, batchSize);
		String[] commandsLin = {ExecutableRunner.findExecutable("pymol"), "-cqr", "PyDescriptor.py"};

		exeRunner.executeBinary(commandsLin, OUTPUFILE, configuration.getTimeoutInSeconds());
		return readResults(OUTPUFILE);
	}

	protected DataTable readResults(String filename) throws IOException
	{
		DataTable dtResults=getResults();

		try{
			BufferedReader in=getAliasedBufferedReader(filename);

			String name =  in.readLine();
			String value = in.readLine();

			in.close();

			String columns[] = name.split(",");
			String values[] = value.split(",");

			if(columns.length != values.length) throw new IOException("Not equal number of descriptors " + values.length + " and descriptor names " + columns.length);

			dtResults.addRow();

			for(int i=1;i<columns.length;i++)
				dtResults.setValue(columns[i], Double.valueOf(values[i]));
		}catch(Exception e){
			dtResults.getCurrentRow().setError(e.getMessage());
		}

		return 	dtResults;	
	}

	@Override
	int getBatchSize(){
		return 1;
	}

	public PyDescriptorsServer() 
	{
		supportedTaskType = DescriptorsConfiguration.PyDescriptor;
		repostSize = 50;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}


}
