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
import qspr.metaserver.configurations.DescriptorsMAP4Configuration;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

public class MAP4DescriptorsServer extends DescriptorsAbstractExecutableServer{

	@Override
	protected
	DataTable calculateDescriptors(DataTable dtMolecules,
			DescriptorsAbstractConfiguration configuration, int start,
			int size) throws Exception {


		DescriptorsMAP4Configuration conf = (DescriptorsMAP4Configuration) configuration;

		saveMolecules(dtMolecules, datain , QSPRConstants.SMILES_FORMAT, start, size);

		String[] commands = new String[] {QSPRConstants.RDKITPYTHON, getExeFile(), "-i", datain, "-o", dataout ,  "--dimensions", ""+ conf.dimensions,"--radius",  ""+conf.radius};

		runPython(commands, dataout, null, dtMolecules.getRowsSize() > 1000 ? dtMolecules.getRowsSize()/100 : dtMolecules.getRowsSize() > 50 ? dtMolecules.getRowsSize() : 50);
		return readResults(dataout,dtMolecules,start,size,configuration);
	}

	@Override
	int getBatchSize(){
		return 100;
	}

	DataTable readResults(String resultFile,DataTable dtMolecules,int start, int size, DescriptorsAbstractConfiguration configuration)throws IOException {

		DescriptorsMAP4Configuration conf = (DescriptorsMAP4Configuration) configuration;

		BufferedReader input = null;
		DataTable res=getResults();
		if(res.getColumnsSize() == 0)
			for(int i=0;i<conf.dimensions;i++)
				res.addColumn("map4_"+i);

		try {
			input = getAliasedBufferedReader(resultFile);

			for(int i=0;i<size;i++) {
				if(dtMolecules.getRow(i+start).isError()) {
					res.addRow().setError(dtMolecules.getRow(i+start).detailedStatus);
					continue;
				}

				String dataResults=input.readLine();
				res.addRow();
				// now we read first line (1 molecule) with results
				String [] initial = dataResults.split("\\s+");
				String [] descriptors = initial[2].split(";");
				if(descriptors.length != conf.dimensions)
					res.getCurrentRow().setError("MAP4 failed");

				for(int l = 0; l < descriptors.length; l++)
					res.setValue(l, Math.log10(Double.parseDouble(descriptors[l])));
			}
		}catch(IOException e) {
			throw(e);
		}
		finally {
			if(input!=null)input.close();
		}

		return res;
	}

	public MAP4DescriptorsServer()
	{
		supportedTaskType = DescriptorsConfiguration.MAP4;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}


}
