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

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsSIGMAConfiguration;
import qspr.workflow.datatypes.DataTable;

public class SIGMAServer extends MOPAC7Server{

	final static String output = "sigma.txt";
	String workingPython;

	public SIGMAServer()
	{
		supportedTaskType = DescriptorsConfiguration.SIGMA;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		repostSize = 10;
		DELIMITER ="\\s+";
	}

	@Override
	String getKeywords(){
		return "EPS=999.0 COSWRT RSOLV=1.4 XYZ PM6 EXTERNAL=POA1.rm1 NSPA=362VDW(H=1.416:C=2.006:N=1.829:O=1.7936:F=1.7346:S=2.124:P=2.124:Cl=2.065:Br=2.183:I=2.3364) GNORM=0.1 RELSCF=0.1" + " THREADS="+THREADS;
	}

	@Override
	String getOutputFile(){
		return FILENAME+".out";
	}

	@Override
	public void parseMopacOutput(String mopacFileName, DataTable dtResult, DescriptorsAbstractConfiguration conf) throws Exception{
		DescriptorsSIGMAConfiguration config = (DescriptorsSIGMAConfiguration) conf;


		boolean first = true;

		for(int nn = 0; nn<DescriptorsSIGMAConfiguration.OPTION_WIDTHS.length; nn++) {

			double width = DescriptorsSIGMAConfiguration.OPTION_WIDTHS[nn];

			if( (config.all == null || !config.all) && width != config.getWidth())continue; // only one width is requested

			String[] commands = {workingPython,"process-cosmo.py","" + width};
			if(workingPython == null)
				workingPython = exeRunner.findWorkingPython(commands, output, 0, null, 10);
			else
				executeBinary(commands, output, 10);

			BufferedReader outputReader = null;
			if(first)dtResult.addRow(); // we add only one time; after that we reuse previous row

			try {
				outputReader = getAliasedBufferedReader(output);

				String line;
				while ((line = outputReader.readLine()) != null)
				{
					line = line.trim();
					if(line.startsWith("#") || line.length() == 0)continue;

					String[] n = line.split(DELIMITER);
					dtResult.setValue("sigma"+nn+"_"+n[0],n[1]);
				}

			}catch(Exception e){
				dtResult.getCurrentRow().setError(e.getMessage());
			}
			finally {
				if(outputReader !=null)outputReader.close();
			}

			first = false;
			System.out.println("**** "+ dtResult.getColumnsSize());
		}

		System.out.println(dtResult.getColumnsSize());

	}


}
