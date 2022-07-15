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

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsGSFragConfiguration;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;

public class GSFragServer extends DescriptorsAbstractExecutableServer {
	static String gsfraglBin;
	static final String gsout="descriptors.txt";

	@Override
	public void setParam(String name, String value) {
		super.setParam(name, value);
		if ("EXE1".equals(name))
			gsfraglBin = value;
	}

	@Override
	int getBatchSize()
	{
		return 1;
	}

	@Override
	public void setStatus(String status) 
	{
		if (status == null || status.contains("Molecule 1"))return;
		super.setStatus(status);
	}

	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration,int start,int batchsize) throws Exception
	{
		if (receivedConfiguration == null || !(receivedConfiguration instanceof DescriptorsGSFragConfiguration))
			throw new Exception("Invalid configuration passed, should be instance of GSFragConfiguration");
		DescriptorsGSFragConfiguration configuration = (DescriptorsGSFragConfiguration) receivedConfiguration;

		saveMolecules(dtMolecules, datain, QSPRConstants.SDF, start, batchsize);

		DataTable dtres = new DataTable(true);

		String[] commands = null;
		if(configuration.isGSFRAG()){
			commands=new String[]{getExeFile(), datain};
			executeBinary(commands, gsout, receivedConfiguration.getTimeoutInSeconds());
			readStandardOutputResults(dtres,gsout,true);
		}

		if(batchsize == 1 && dtres.getRowsSize()  == 1 && dtres.getRow(0).isError()){ 
			//  If previous calculations failed - row marked as failed and second part of calculations is skipped to speed up calculations ...
			DataTable dtResult = getResults(); // previous result table
			dtResult.addRow(dtres.getRow(0));
			return dtResult;
		}

		if(configuration.isGSFRAGL()){
			DataTable dtres1 = new DataTable(true);
			commands=new String[]{getAliasedFileName(gsfraglBin), datain};
			executeBinary(commands, gsout, receivedConfiguration.getTimeoutInSeconds());
			readStandardOutputResults(dtres1,gsout,true);
			if(configuration.isGSFRAG())dtres.mergeColumnsWith(dtres1);
			else
				dtres = dtres1;
		}

		DataTable dtResult = getResults(); // previous result table

		// now we will add all descriptors to previous result
		// one by one row and considering that new column may have appear
		dtres.reset();
		while(dtres.hasMoreRows()){
			dtres.nextRow(); 
			dtResult.addRow();
			if ("error".equals(dtres.getCurrentRow().status))
				dtResult.getCurrentRow().setError(dtres.getCurrentRow().detailedStatus);
			else
				for(int i=0;i<dtres.getColumnsSize();i++)
					dtResult.setValue(dtres.getColumn(i),dtres.getValue(i));
		}

		return dtResult;
	}

	public GSFragServer() 
	{
		repostSize = 200;
		supportedTaskType = DescriptorsConfiguration.GSFrag;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}
}
