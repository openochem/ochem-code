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

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Map;

import qspr.metaserver.configurations.LinearConfiguration;
import qspr.metaserver.configurations.MLRAConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.utils.FileUtils;

public class MLRAServer extends LinearAbstractServer
{
	static int MAXSTEPS = 33;

	@Override
	protected Map<Integer, Double> calculateCoefficients(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, LinearConfiguration configuration) throws Exception{

		BufferedWriter writer = getAliasedBufferedWriter(DATAFILE);
		writeSimpleIgorFormatData(dtDescriptors, dtExpValues, writer);
		writer.close();
		MLRAConfiguration config = (MLRAConfiguration)configuration;
		try{
			runBinary(config,dtExpValues.getDataSize());
		}catch(IOException e){
			throw new CriticalException("MLRA Could not find optimal parameters. Either dataset is too small or descriptors are too much correlated.");
		};
		setStatus("MLRA has finished - reading results.");
		return getCoefficients(MODEL,false);
	}

	private void runBinary(MLRAConfiguration configuration, int size) throws IOException {
		int nonzerolimit= size/2 > 1000 ? 1000 : size/2;
		double correlation=1;

		for(int steps=1, nonzero=1;steps<=MAXSTEPS;steps++)
			try{
				createCfgFile(configuration, nonzero, correlation);
				setStatus("optimizing mlra step "+steps+" ouf of "+MAXSTEPS+" NONZERO "+nonzero+" correlation "+correlation);

				nonzero *= 2;
				if(nonzero > nonzerolimit)nonzero=nonzerolimit;
				correlation *= 0.95;

				executeBinary();
				if(new File(getAliasedFileName(OCHEM)).exists())return;

			}catch(Exception e)
		{
				out.println("step: "+steps+" error: "+e.getMessage());
		}

		try{
			createCfgFile(configuration, nonzerolimit, 0.01);
			executeBinary();
			if(new File(getAliasedFileName(OCHEM)).exists())return;

		}catch(Exception e)
		{
			out.println("failed totally");
		}


		throw new IOException("MLRA Could not find optimal parameters in runBinary.");
	}

	/**
	 * MLRA cfg file template FILE=data NONZERO=1 CORRELATION=1 MODELS=1 NORM=1
	 * NAMES=1 ANALYSIS=1
	 * 
	 * ALPHA=0.05 STOP END
	 */
	private void createCfgFile(MLRAConfiguration configuration,int nonzero,double correlation)
			throws IOException {
		String parameters = "FILE=data\nNAMES=1\nNORM=1\nMODELS=1\nANALYSIS=1\n";
		parameters+="\nNONZERO="+nonzero
				+"\nCORRELATION="+correlation
				+ "\nALPHA="+ configuration.alpha 
				+ "\nNVAR="+ configuration.nvariables 
				+ "\nSTOP\nEND";

		FileUtils.saveStringToFile(parameters, getAliasedFileName(CFG));
	}

	public MLRAServer() {
		supportedTaskType = QSPRConstants.MLRA;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setInputFlowGroup(2);
		setOutputFlowGroup(0);
	}


}
