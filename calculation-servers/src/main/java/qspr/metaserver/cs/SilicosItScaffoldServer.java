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
import java.io.IOException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsSilicosItScaffoldConfiguration;
import qspr.metaserver.util.CanonicalSmiles;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;


public class SilicosItScaffoldServer extends DescriptorsAbstractExecutableServer 
{

	private static transient final Logger logger = LogManager.getLogger(SilicosItScaffoldServer.class);

	public SilicosItScaffoldServer() 
	{
		supportedTaskType = DescriptorsConfiguration.SilicosItScaffold; // Should be a unique string - "task type"
		repostSize = 1000; // The dataset will be automatically broken to pieces of this size and distributed for parallel calculations 

		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}

	@Override
	public int getBatchSize() 
	{
		return 1; // The dataset will be automatically broken to pieces of this size within one calculation server and calculated consequently 
	}


	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration configuration, int start, int batchSize) throws Exception 
	{
		DescriptorsSilicosItScaffoldConfiguration conf;
		if (configuration == null) {
			//TODO default configuration or throw an exception?
			conf = new DescriptorsSilicosItScaffoldConfiguration();
		} else {
			conf = (DescriptorsSilicosItScaffoldConfiguration) configuration;
		}

		int size = (dtMolecules.getRowsSize() > start + batchSize) ? batchSize : dtMolecules.getRowsSize()-start;

		return calculateStripitScaffolds(dtMolecules, conf, start, size); // A sample of the calculation server with pure java
	}


	private DataTable calculateStripitScaffolds(DataTable dtMolecules, DescriptorsSilicosItScaffoldConfiguration configuration, int start, int size) throws Exception {
		logger.info("Calculating with external binary");
		DataTable dtResults = getResults();

		if(size != 1) throw new IOException("Can process only one molecule per batch.");
		dtResults.addRow(); // only one row is added (and it is always added)

		// We prepare a file with molecules for an external tool
		final String inputFile = "data.sdf";
		final String outputFile = "result";
		final String configurationFile = "configuration";

		saveMolecules(dtMolecules, inputFile, QSPRConstants.SDF, start, size);
		BufferedWriter bw = getAliasedBufferedWriter(configurationFile);
		configuration.writeScaffoldParameters(bw);
		bw.close();

		// Form a command line for the external tool
		String[] commands = {getExeFile(), "--input", getAliasedFileName(inputFile), "--scaffolds", configurationFile, "--output", outputFile, "--noLog"};

		// Execute a tool.... we handle timeouts and other stuff like that
		executeBinary(commands, outputFile,10*size);

		//Read the output
		BufferedReader reader = null;

		try {
			reader = getAliasedBufferedReader(outputFile); //XXX use of AliasedBufferReader, is it correct?
			String s = reader.readLine(); //Headers
			String[] headerPieces = s.split("\\s+");
			s = reader.readLine(); // Scaffolds

			String[] bodyPieces = s.split("\\s+");

			if (bodyPieces.length != headerPieces.length) {
				dtResults.getCurrentRow().setError("Number of descriptors does not match number of headers: " + bodyPieces.length + " != " + headerPieces.length);
				return dtResults;
			}

			if (bodyPieces[2].charAt(0) == '-') {
				dtResults.setValue(QSPRConstants.SMILES_FORMAT+":"+"C", 1.); // If there is no ring, Strip-it writes '-' to file
				return dtResults;
			}

			for (int i = 2; i<headerPieces.length; i++) {
				dtResults.setValue(QSPRConstants.SMILES_FORMAT+":"+CanonicalSmiles.canonize(bodyPieces[i]), 1.); // Strip-it returns canonical smiles	
			}
		}catch(Exception e) {
			throw e;
		}finally {
			reader.close();
		}

		return dtResults;

	}

}
