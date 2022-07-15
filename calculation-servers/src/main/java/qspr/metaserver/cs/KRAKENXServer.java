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


import java.io.IOException;

import com.eadmet.utils.FileUtils;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

public class KRAKENXServer extends MOPAC7Server{

	final static String output = "mols_desc.txt";

	//String ignore[] = {"MOPAC_POINTGROUP","MOPAC_POL_ALPHA","MOPAC_POL_BETA","MOPAC_POL_GAMMA","CPSA_RPCG"};

	public KRAKENXServer()
	{
		supportedTaskType = DescriptorsConfiguration.KRAKENX;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		repostSize = 20;
		DELIMITER ="\\s+";
	}

	@Override
	String getKeywords(){
		return "PM6 PRECISE MMOK XYZ ENPART STATIC SUPER BONDS LARGE PRTXYZ DISP ALLVEC THREADS="+THREADS;
	}

	@Override
	String getOutputFile() throws IOException{
		FileUtils.saveStringToFile("FOR005.out",getAliasedFileName("mopac.txt"));
		FileUtils.saveStringToFile("FOR005.sdf",getAliasedFileName("mols.txt"));
		return FILENAME+".out";
	}

	@Override
	public void parseMopacOutput(String mopacFileName, DataTable dtResult, DescriptorsAbstractConfiguration conf) throws IOException, InterruptedException{
		String[] commands = { "cd",getAliasedFileName(".")+";",javaHome+"/bin/java", "-Djava.awt.headless=true", "-Duser.language=en", "-Duser.region=US",
				"-cp",".:"+javaClassPath+":KrakenX.jar", // we need CDK2 ...
				"krakenx.KrakenX","pars.txt"};
		executeBinaryBash(commands, output, TIMEOUTBABEL);
		DataTable test = new DataTable(true);
		
		String ou = FileUtils.getFileAsString(getAliasedFileName(output));
		if(!ou.contains("molecule"))FileUtils.saveStringToFile("molecule "+ou, getAliasedFileName(output)); // fixing recent bug in Krakenx, removal of the name of first column
		
		readStandardOutputResults(test, output,false);
		readStandardOutputResults(dtResult, output,false);

		AbstractDataRow r = dtResult.getCurrentRow();

		if(r.isError()) {
			r.setError(r.detailedStatus+ " in data input");
			return;
		}

		if(test.getColumnsSize() != 124) {
			r.setError("number of columns != 124 (expected number)");
			return;
		}

		int n = 0;
		for(String s:dtResult.columns) 
			if(((Double)dtResult.getValue(s)).isInfinite() || ((Double)dtResult.getValue(s)).isNaN()) {
				dtResult.setValue(s,0);
				n++;
			}

		if(n == dtResult.getColumnsSize()) {
			r.setError("all NaN in data input");
		}

	}


}
