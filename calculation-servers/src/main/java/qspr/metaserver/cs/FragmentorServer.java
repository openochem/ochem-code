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
import java.io.File;
import java.io.IOException;
import java.util.Vector;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsFragmentorConfiguration;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;

public class FragmentorServer extends DescriptorsAbstractExecutableServer
{
	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration,
			int start,int size) throws Exception
	{
		DescriptorsFragmentorConfiguration configuration=(DescriptorsFragmentorConfiguration)receivedConfiguration;
		String file = (new File(getExeFile())).exists()?getExeFile():getAliasedFileName("isida-fragmentor-linux");
		String commands[]={file,
				"-i",datain,
				"-o",dataout,
				"-t",""+configuration.fragmentType,
				"-l",""+configuration.minFragmentLength,
				"-u",""+configuration.maxFragmentLength,
				"-t","6",
				"-t","7",
				"-f","SVM"};

		saveMolecules(dtMolecules, datain,QSPRConstants.SDFAROM_BASIC_WITHH,start,size);
		executeBinary(commands, dataout+".svm", 15);
		return readResults(dataout+".svm",dtMolecules.getRowsSize());
	}

	@Override
	int getBatchSize(){
		return 1000;
	}

	public FragmentorServer()
	{
		supportedTaskType = DescriptorsConfiguration.FRAGMENTS;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}

	DataTable readResults(String resultFile,int totalMolSize)throws IOException {
		Vector<String> descriptorNames = new Vector<String>();

		// First we have to read names of descriptors
		String temp = null;
		BufferedReader bis = getAliasedBufferedReader(dataout+".hdr");
		while((temp = bis.readLine()) != null)
		{
			String[] descriptors = temp.trim().split("\\s+");
			descriptorNames.add(descriptors[1]);
		}
		bis.close();

		BufferedReader input = getAliasedBufferedReader(resultFile);
		DataTable res=getResults();
		String dataResults=null;

		while((dataResults=input.readLine())!=null){
			res.addRow();
			// now we read first line (1 molecule) with results
			String [] descriptors = dataResults.split("\n")[0].trim().split("\\s+");
			for(int l = 1; l < descriptors.length; l++)
			{
				String []desc = descriptors[l].split(":");
				double value = (double)Integer.parseInt(desc[1]);
				String name=descriptorNames.get(Integer.parseInt(desc[0])-1);
				if(name==null)throw new IOException("Descriptor was not found by index");
				// we are limiting number of descriptors x name to the maximum allowed limit per table
				//				if(totalMolSize*res.getColumnsSize()<QSPRConstants.DESCRIPTORS_X_MOLECULES
				//						|| 	res.containsColumn(name))
				res.setValue(name, value);
			}
		}
		input.close();
		return res;
	}

}
