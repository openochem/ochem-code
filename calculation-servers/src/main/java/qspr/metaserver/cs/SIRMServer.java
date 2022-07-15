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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;

import qspr.dao.ChemDAO;
import qspr.dao.Various;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsSIRMSConfiguration;
import qspr.metaserver.protocol.Task.DebugLevel;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

public class SIRMServer extends DescriptorsAbstractExecutableServer{

	static String workingPython;

	@Override
	protected
	DataTable calculateDescriptors(DataTable dtMolecules,
			DescriptorsAbstractConfiguration configuration, int start,
			int batchSize) throws Exception {

		DescriptorsSIRMSConfiguration conf = (DescriptorsSIRMSConfiguration)configuration;

		if(conf.minFragment <= 0 ) throw new IOException("Incorrect configuratin parameters Minimum Fragment " + conf.minFragment + " <= 0");

		if(conf.minFragment > conf.maxFragment) throw new IOException("Incorrect configuratin parameters: min fragmnet = " + conf.minFragment + " < max fragment = " + conf.maxFragment);

		DataTable mol = new DataTable();
		mol.addColumn(QSPRConstants.SDF_COLUMN);

		setStatus("starting adding properties");

		for(int i=start; i<start+batchSize;i++)try{
			String sdf = (String)dtMolecules.getValue(i, 0);
			sdf = Various.molecule.convertToFormat(sdf, QSPRConstants.SDFAROM_BASIC_WITHH);

			sdf = sdf.replace("$$$$\n", "");
			
			for(DescriptorsSIRMSConfiguration.Labeling l:conf.getUniqueAndSorted())
				switch(l){

				case HB:   sdf = Various.molecule.addAtomProperty(ChemDAO.Properties.HB, sdf); break; // in fact it is not used since not defined in setup.txt
				case LOGP: sdf = Various.molecule.addAtomProperty(ChemDAO.Properties.LOGP, sdf); break;
				case REFRACTIVITY: sdf = Various.molecule.addAtomProperty(ChemDAO.Properties.REFRACTIVITY, sdf); break;
				case CHARGE: sdf = Various.molecule.addAtomProperty(ChemDAO.Properties.CHARGE, sdf); break;
				case elm:
				case none: break;
				}
			mol.addRow();
			mol.setValue(sdf);	
		}catch(Exception e) {
			setStatus("exception " + e.getMessage());
			if(batchSize > 1) throw e;
			DataTable dtResults = getResults();
			if(debug != DebugLevel.NONE)
				e.printStackTrace();
			dtResults.addRow().setError("SIRMS: failed to process this molecule");
			return dtResults;
		}

		setStatus("saving molecules");

		saveMolecules(mol, datain + ".sdf", conf.ignoreH() ?QSPRConstants.SDFNOH : QSPRConstants.SDFH, 0, batchSize);

		if(conf.minFragment > conf.maxFragment) throw new IOException("Min value " + conf.minFragment + "  > " + conf.maxFragment);

		String[] commands = new String[] { workingPython, getExeFile(), "-i", datain+".sdf", "-o", dataout,
				"--min_atoms",""+conf.minFragment,"--max_atoms",""+conf.maxFragment,"--output_format","svm"
		};

		if(conf.asMixtures())commands = OCHEMUtils.append(commands, "-q");

		if(conf.descriptorTypes.size()>0){
			if(conf.ignoreH())commands = OCHEMUtils.append(commands, "-a");
			for(DescriptorsSIRMSConfiguration.Labeling l:conf.getUniqueAndSorted())
				commands = OCHEMUtils.append(commands, ""+l);
		}

		if(workingPython == null)
			workingPython = exeRunner.findWorkingPython(commands, dataout, 0, null, dtMolecules.getRowsSize() > 1000 ? dtMolecules.getRowsSize()/100 : dtMolecules.getRowsSize());
		else
			executeBinary(commands, dataout);

		return readResults(dataout);
	}

	@Override
	int getBatchSize(){
		return 10;
	}

	DataTable readResults(String filename) throws IOException
	{
		DataTable dtResults = getResults();

		//map new columns to old ones
		BufferedReader in=getAliasedBufferedReader(filename+".colnames");
		if(in == null)throw new UserFriendlyException("No file was produced: " + filename+".colnames");
		Map<Integer,Integer> newColToResCol = new HashMap<Integer,Integer>();
		List<String> newcols = new ArrayList<String>();
		String dataResults=null;
		for(int i=0; (dataResults=in.readLine())!=null; i++){
			newcols.add(dataResults);
			int index = dtResults.getColumnIndex(dataResults);
			if(index == -1){
				dtResults.addColumn(dataResults);
				index = dtResults.getColumnIndex(dataResults);
			}
			newColToResCol.put(i, index);
		}
		in.close();

		File f = new File(getAliasedFileName(filename+".colnames")); f.delete();

		in=getAliasedBufferedReader(filename);
		if(in == null)throw new UserFriendlyException("No file was produced: " + filename);

		while((dataResults=in.readLine())!=null){

			boolean oneMol=false;

			dtResults.addRow();

			String[] pieces = dataResults.trim().split("\\s+");
			for(String piece:pieces){
				String[] vals = piece.split(":");
				dtResults.setValue(newColToResCol.get(Integer.valueOf(vals[0])), vals[1]);
				oneMol = true;
			}

			if(!oneMol)
				dtResults.getCurrentRow().setError("SIRMS did not calculate any descriptors for this molecule.");
		}			

		in.close();

		return dtResults;
	}

	public SIRMServer()
	{
		repostSize = 1000;
		supportedTaskType = DescriptorsConfiguration.SIRMS;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}

}
