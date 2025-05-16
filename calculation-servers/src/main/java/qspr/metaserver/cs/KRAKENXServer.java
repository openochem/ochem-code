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
import com.eadmet.utils.OSType;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

public class KRAKENXServer extends MOPAC7Server{

	final static String output = "mols_desc.txt";

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

		String java = "/usr/lib/jvm/java-11-openjdk-"+
				(OSType.isAarch64()?"arm64":"amd64")+
				"amd64/bin/java";

		String[] commands = { "cd",getAliasedFileName(".")+";",
				java, "-Djava.awt.headless=true", "-Duser.language=en", "-Duser.region=US",
				"-cp",".:"+javaClassPath+":KrakenX.jar", // we need CDK2 ...
				"krakenx.KrakenX","pars.txt"};
		FileUtils.saveStringToFile("FOR005.out", getAliasedFileName("mopac.txt"));
		executeBinaryBash(commands, output, TIMEOUTBABEL);
		DataTable test = new DataTable(true); 
		String out = FileUtils.getFileAsString(getAliasedFileName(output));
		if(!out.contains("molecule")) { // fix for bug with shift of columns
			out = "molecule " + out;
			FileUtils.saveStringToFile(out, getAliasedFileName(output));
		}

		if(out.contains(" NA ")) { // fix for bug with shift of columns
			out=out.replaceAll(" NA ", " NaN ");
			out=out.replaceAll(" NA ", " NaN ");
			FileUtils.saveStringToFile(out, getAliasedFileName(output));
		}

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
