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
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.LineNumberReader;

import com.eadmet.utils.OSType;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

public class CDDDDescriptorsServer extends DescriptorsAbstractExecutableServer
{

	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration,
			int start,int size) throws Exception
	{

		DataTable res=getResults();

		try {
			int mol = prepareData(dtMolecules, start, size);
			if(mol > 0) {
				String commands[]={"export","PATH=$PATH:~/.local/bin",";","cddd-onnx","--input",getAliasedFileName(datain+".smi"),"--output",getAliasedFileName(dataout)};
				runPython(commands, dataout, CONDA.RDKIT, 50);
				return readResults(dataout,dtMolecules,start,size);
			}
		}catch(Throwable e) {
			if(size > 1)throw e;
		}

		res.addRow().setError("CDDD calculations failed");
		return res;
	}

	int prepareData(DataTable dtMolecules, int start, int size) throws Exception{

		saveMolecules(dtMolecules,QSPRConstants.SMILES_FORMAT,QSPRConstants.SMILESNOAROM, start, size);
		String filePython = (new File(QSPRConstants.PYTHON36)).exists()?QSPRConstants.RDKITPYTHON:"python3";
		String[] commands = new String[] {OSType.isMac()?"python3":filePython, getExeFile(), "--infile", QSPRConstants.SMILES_FORMAT, "--outfile",
				dataout,"--augment","1","--isomeric","True"};
		exeRunner.findWorkingPython(commands,dataout,0,CONDA.RDKIT, size);

		LineNumberReader br =  getAliasedLineNumberReader(dataout);
		BufferedWriter writer  =  getAliasedBufferedWriter(datain+".smi");

		int count = 0;
		String line;
		for(int row = start; (line = br.readLine()) != null;row++) {
			String smiles[] = line.split(",");
			if(smiles.length != 2)continue;
			if(smiles[1].indexOf(QSPRConstants.ERROR) != -1) {
				dtMolecules.getRow(row).setError("RDKiT failed to convert this molecule :" + smiles[0]);
				continue;
			}
			writer.append(smiles[1]+"\n");
			count++;
		}

		writer.append("O=CCCNCO\n"); // one more -- just to always have correct output
		br.close();
		writer.close();

		return count;
	}


	@Override
	int getBatchSize(){
		return 100;
	}

	public CDDDDescriptorsServer()
	{
		supportedTaskType = DescriptorsConfiguration.CDDD;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		repostSize = 1000;
	}

	DataTable readResults(String resultFile,DataTable dtMolecules,int start, int size)throws IOException {

		BufferedReader input = null;
		DataTable res=getResults();
		if(res.getColumnsSize() == 0)
			for(int i=1;i<=512;i++)
				res.addColumn("cddd_"+i);

		try {
			input = getAliasedBufferedReader(resultFile);
			String dataResults=input.readLine();

			for(int i=0;i<size;i++) {
				if(dtMolecules.getRow(i+start).isError()) {
					res.addRow().setError(dtMolecules.getRow(i+start).detailedStatus);
					continue;
				}

				dataResults=input.readLine();
				res.addRow();
				// now we read first line (1 molecule) with results
				String [] descriptors = dataResults.split(",");
				if(descriptors.length != 514 || descriptors[2].length()==0) {
					res.getCurrentRow().setError("CDDD failed");
					continue;
				}

				for(int l = 2; l < descriptors.length; l++)
					res.setValue(l-2, descriptors[l]);
			}
		}catch(IOException e) {
			throw(e);
		}
		finally {
			if(input!=null)input.close();
		}

		return res;
	}

}
