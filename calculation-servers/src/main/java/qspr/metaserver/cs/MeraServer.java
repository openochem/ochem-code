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
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;

public class MeraServer extends DescriptorsAbstractExecutableServer
{

	/**
	 *  TimeOut for calculation of one molecules
	 * @return
	 */
	int getTimeOut(int size){
		return size==1?10:size*3;
	}

	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration, int start, int size) throws Exception
	{
		String[] commands = {getExeFile(), datain};
		saveMolecules(dtMolecules, datain, QSPRConstants.SDF, start, size);
		try {
			String res=datain+"."+supportedTaskType.toLowerCase();
			executeBinary(commands, res, getTimeOut(size));
			return readResults(res);
		}catch(Exception e) {
			if(size != 1)throw e;
			DataTable dtResults=getResults();
			dtResults.addRow();
			dtResults.getCurrentRow().setError(supportedTaskType + " did not calculate any descriptors for this molecule.");
			return dtResults;
		}
	}

	@Override
	int getBatchSize(){
		return 100;
	}

	DataTable readResults(String filename) throws IOException
	{
		DataTable dtResults=getResults();
		BufferedReader in=getAliasedBufferedReader(filename);
		String dataResults=null;

		boolean oneMol=false;

		for(int l=0;(dataResults=in.readLine())!=null;l++){

			String[] pieces = dataResults.trim().split("\\s+");

			if (dtResults.getColumnsSize() == 0)
				for(String descriptor: pieces)
					dtResults.addColumn(descriptor);

			if(l==0) continue; // Skip header

			dtResults.addRow();
			oneMol=true;

			if (pieces.length != dtResults.getColumnsSize())
			{
				dtResults.getCurrentRow().setError(supportedTaskType + " could not process this compound - number of descriptors did not match the expected number");
				continue;
			}

			for (int i=0; i<pieces.length; i++)
				dtResults.setValue(i, pieces[i].trim());
		}
		in.close();

		if(!oneMol){
			dtResults.addRow();
			dtResults.getCurrentRow().setError(supportedTaskType + " did not calculate any descriptors for this molecule.");
		}

		return 	dtResults;	
	}
	

	public MeraServer()
	{
		supportedTaskType = DescriptorsConfiguration.Mera;
		repostSize = 10000;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}
}
