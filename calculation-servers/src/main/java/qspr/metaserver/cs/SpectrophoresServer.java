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
import qspr.metaserver.configurations.DescriptorsSpectrophoresConfiguration;
import qspr.metaserver.util.ExecutableRunner;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;

public class SpectrophoresServer extends DescriptorsAbstractExecutableServer {

	String SpectrophoresNames[]={"SpectrophoresPartial","SpectrophoresLipophilicity","SpectrophoresShape","SpectrophoresElectrophilicity"};

	@Override
	protected
	DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration, int start, int batchSize)
	throws Exception {

		DescriptorsSpectrophoresConfiguration configuration = (DescriptorsSpectrophoresConfiguration) receivedConfiguration;

		saveMolecules(dtMolecules, datain+".sdf",QSPRConstants.SDF,start,batchSize);
		String[] commands = new String[] { ExecutableRunner.findExecutable("obspectrophore"), "-i", getAliasedFileName(datain+".sdf"), 
				"-a",""+configuration.getAccuracy(),
				"-s",configuration.getCageTypeString(),
				"-r",""+configuration.getResolution()};

		executeBinary(commands, stdout);

		int cagetype=configuration.getCageType();
		
		return readResults(stdout,cagetype==DescriptorsSpectrophoresConfiguration.cageTypeNo?48
				:cagetype==DescriptorsSpectrophoresConfiguration.cageTypeUnique?72
						:cagetype==DescriptorsSpectrophoresConfiguration.cageTypeMirror?72
								:cagetype==DescriptorsSpectrophoresConfiguration.cageTypeAll?144:0
					);
	}

	@Override
	int getBatchSize(){
		return 100;
	}
	
	DataTable readResults(String result,int n) throws IOException{
		DataTable dtResults=getResults();

		boolean first=true;
		
		// If there are no names yet, we generate them
		if(dtResults.getColumnsSize()==0){
			if(n==48)
			for(int i=0;i<48;i++)
				dtResults.addColumn(SpectrophoresNames[i/12]+"_"+(1+(i%12)));
			else
				for(int i=0;i<n;i++)
					dtResults.addColumn("Spectrophores_"+i);
		}
		
		BufferedReader reader=getAliasedBufferedReader(result);
		
		String line=null;
				
		while((line=reader.readLine())!=null){
			String values[]=line.split("\\s+");
			if(values.length<48)continue; // to skip lines without any data		
			dtResults.addRow();
			first=false;
			for(int i=1;i<values.length;i++){
				if(values[i].toLowerCase().equals("nan")){
					dtResults.getCurrentRow().setError("Spectrophores failed for "+dtResults.getColumn(i-1));
					break;
				}
				dtResults.setValue(i-1, Double.parseDouble(values[i]));
			}
		}
		reader.close();
		
		if(first){
			dtResults.addRow();
			dtResults.getCurrentRow().setError("Spectrophores failed to calculate any descriptors for this molecule");
		}
			
		return dtResults;
	}
	
	public SpectrophoresServer()
	{
		supportedTaskType = "Spectrophores";
		repostSize = 5000;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}

}
