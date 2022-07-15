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
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsRDKITConfiguration;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;

public class RDKITDescriptorsServer extends DescriptorsAbstractExecutableServer{

	protected static final String CFG = "config.cfg";
	protected static final String INPUT = "input.sdf";
	protected static final String OUTPUT = "descriptors.txt";
	String python = null;

	@Override
	protected DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration configuration, int start, int batchSize)
			throws Exception {
		DescriptorsRDKITConfiguration conf = (DescriptorsRDKITConfiguration) configuration;

		createCfg(conf);
		System.out.print("running "+ conf);
		saveMolecules(dtMolecules, INPUT, QSPRConstants.SDF, start, batchSize);
		String[] commandsLin = {python, getExeFile(), getAliasedFileName(CFG)};
		runPython(commandsLin, OUTPUT, CONDA.RDKIT, batchSize>100?batchSize/10:10);
		return readResults(OUTPUT);
	}


	String strim(String value){
		return value.replaceAll("['{},:\"]", "").toUpperCase();
	}

	protected DataTable readResults(String filename, DataTable dtResults) throws IOException
	{
		BufferedReader in=getAliasedBufferedReader(filename);

		Set<String> ignore = new HashSet<String>(Arrays.asList("cansmi", "cansmiNS", "formula", "InChI","InChIKey", 
				"L5", "s", "smarts"));

		boolean oneMol=false;
		String dataResults=null;
		while((dataResults=in.readLine())!=null){

			if(!dataResults.startsWith("{"))continue;

			String[] pieces = dataResults.trim().split("\\s+"); 

			if(pieces.length<=2)continue; // there are no descriptors, skip

			dtResults.addRow();
			oneMol=true;

			for(int i=0;i<pieces.length;i+=2){
				String name=strim(pieces[i]);
				if(ignore.contains(name))continue;
				addValueIgnoreNaN(dtResults,name,strim(pieces[i+1]));
			}
		}
		in.close();

		if(!oneMol){
			dtResults.addRow();
			dtResults.getCurrentRow().setError("RDKit did not calculate any descriptors for this molecule.");
		}

		return 	dtResults;	
	}

	void addValueIgnoreNaN(DataTable dtResults, String name, String value) {
		if(value.equalsIgnoreCase("NAN")){
			if(dtResults.getColumnIndex(name)==-1)return; // this descriptor has only NaN values, ignore it
			value ="0";
		}
		dtResults.setValue(name,value);
	}

	protected boolean createCfg(DescriptorsRDKITConfiguration configuration) throws IOException {

		BufferedWriter writer = getAliasedBufferedWriter(CFG);
		writer.write("[Task]\n");
		writer.write("\ninput_file = " + INPUT);
		writer.write("\noutput_file = " + OUTPUT);

		int additional = 0;

		for (int i = 1; i <= configuration.getBlocks(); i++)
			if(configuration.bit(i-1))
				switch(i) {
				case 1: writer.write("\nscalars = True");break;
				case 2: writer.write("\nscalars_secondary =  True");break;
				case 3: writer.write("\nautocorr2d = True");break;
				case 4: writer.write("\nautocorr3d = True");break;
				case (DescriptorsRDKITConfiguration.TOPOLOGICAL +1): 
					writer.write("\ntopological = True");
				writer.write("\ntopological_nbits = " + configuration.TOPOLOGICAL_NBITS + "\n");
				break;
				case 6: writer.write("\ngetaway = True");break;
				case 7: writer.write("\nmorse = True");break;			
				case 8: writer.write("\nrdf = True");break;

				case (DescriptorsRDKITConfiguration.WHIM +1): //RDKitnewConfiguration.WHIM 
					writer.write("\nwhim = True");
				writer.write("\nwhim_thresh = " + configuration.WHIM_THRESHOLD + "\n");
				break;
				case (DescriptorsRDKITConfiguration.MORGAN +1): //RDKitnewConfiguration.MORGAN
					writer.write("\nmorgan = True");
				writer.write("\nmorgan_nbits = " + configuration.MORGAN_NBITS);
				writer.write("\nmorgan_radius = " + configuration.MORGAN_RADIUS);
				if(configuration.MORGAN_FCFP != null)writer.write("\nmorgan_fcfp = " + configuration.MORGAN_FCFP + "\n");
				if(configuration.MORGAN_COUNTS != null)writer.write("\nmorgan_counts = " + configuration.MORGAN_COUNTS);
				break;
				case 11: 
					writer.write("\nmaccs = True");
					break;

				case 12: 
					writer.write("\natom_pairs = True");
					break;
				case 13: 
					writer.write("\nbtf = True");
					break;
				case 14: 
					writer.write("\nbpf = True");
					break;
				case 15: 
					writer.write("\ntorsions = True");
					break;
				case 16: writer.write("\nsascore = True");break;
				case (DescriptorsRDKITConfiguration.AVALON +1): //Avalon
					writer.write("\navalon = True");
					writer.write("\navalon_nbits = " + configuration.AVALON_NBITS);
					if(configuration.AVALON_COUNTS != null && configuration.AVALON_COUNTS)writer.write("\navalon_count = True");
					break;
				case 18: writer.write("\nrdkitdef = True");break;
					
				default: throw new IOException("No descriptors are defined for block " + i);
				}

		writer.write("\n");
		writer.close();

		return  additional != configuration.dragonBlocks;
	}


	protected DataTable readResults(String filename) throws IOException
	{
		DataTable dtResults=getResults();
		BufferedReader in=getAliasedBufferedReader(filename);

		String dataResults=null, columns[] = null;

		while((dataResults=in.readLine())!=null){

			String[] pieces = dataResults.trim().split("\\s+"); 

			if(pieces.length == 0 || (pieces.length == 1 && pieces[0].length() == 0) ) continue;

			if(columns == null) {
				columns = pieces;
				continue;
			}

			dtResults.addRow();
			for(int i = 0; i < pieces.length; i++) 
				addValueIgnoreNaN(dtResults,columns[i],pieces[i]);
		}
		in.close();

		return 	dtResults;	
	}

	@Override
	int getBatchSize(){
		return 20;
	}

	public RDKITDescriptorsServer()
	{
		supportedTaskType = DescriptorsConfiguration.RDKIT;
		repostSize = 3000;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}

}
