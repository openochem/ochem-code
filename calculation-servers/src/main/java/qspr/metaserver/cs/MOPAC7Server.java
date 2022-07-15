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
import java.util.Arrays;
import java.util.List;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;

import qspr.dao.ChemDAO.Properties;
import qspr.dao.ChemDAO;
import qspr.dao.ChemInfEngine;
import qspr.dao.Various;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.util.ExecutableRunner;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;

public class MOPAC7Server extends DescriptorsAbstractExecutableServer
{
	/**
	 * MOPAC keywords to use. 0SCF will not work because it only prints data,
	 * then stops. RESTART will not help with this, either. The EF keyword is
	 * not compatible with the 1SCF keyword.
	 */

	public final static String FILENAME = "FOR005";

	static final int TIMEOUTBABEL = 60; // one minute is even probably too much

	static final String HOMO = "HomoEnergy", LUMO ="LumoEnergy",
			KEYWORDS = "XYZ ENPART POLAR PRECISE BONDS LARGE VECTORS=ALLVEC GEO-OK PI MMOK", AM1 ="AM1", PM7 = "PM7";

	static final String ERRORSTANDARD ="Failed and did not produce any output";

	static final String datainsdf = FILENAME + ".sdf";

	protected static final int THREADS = 8;

	public MOPAC7Server()
	{
		supportedTaskType = DescriptorsConfiguration.MOPAC;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
		repostSize = 30;
	}

	String getKeywords(){
		return AM1 + " " + KEYWORDS;
	}

	int getBatchSize()
	{
		return 1;
	}

	@Override
	final public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration,int start,int size) throws IOException 
	{
		setStatus("Running " + supportedTaskType + " " + start);

		DataTable dtResults = getResults();

		for(int mol = start; mol < start + size ; mol++){

			dtMolecules.setCurrentRow(mol);
			setStatus("processing " + mol + " out of " + dtMolecules.getRowsSize() + " molecules");

			String sdf = dtMolecules.getSDF();

			String charge = null;

			try {
				charge = Various.molecule.getAtomProperty(Properties.CHARGE, sdf);
			}catch(Exception e) {
				try{
					ChemDAO cdk = Various.getCheminfImpl(ChemInfEngine.CDK);
					charge = cdk.getAtomProperty(Properties.CHARGE, sdf);
				}catch(Exception ee) {
					throw new CriticalException("CDK upload failed");
				}
			}

			double value =0;
			for(String s:charge.split(";"))
				value += Double.parseDouble(s);

			String add = Math.abs(value) <0.5?"":value >0? "CHARGE=+"+value:"CHARGE="+value;

			boolean canon = false;

			for(int cycle = 0; cycle < 3; cycle++){

				boolean addedOptionalKeyword = false;

				String message = calc(dtMolecules,dtResults,add,canon,receivedConfiguration);

				if(message != null && message.contains("EXCESS NUMBER OF OPTIMIZATION CYCLES"))
				{message = calc(dtMolecules,dtResults,add += " GNORM=1 ",canon,receivedConfiguration); addedOptionalKeyword = true;}

				if(message != null && message.contains("Add keyword RHF and re-run"))
				{message = calc(dtMolecules,dtResults,add += " RHF ",canon,receivedConfiguration); addedOptionalKeyword = true;}

				if(message != null && message.contains("DUE TO PROGRAM BUG"))
				{	
					canon = true;
					message = calc(dtMolecules,dtResults,add,canon,receivedConfiguration); addedOptionalKeyword = true;
				}

				if(message != null && message.contains("COMPONENTS OF ALPHA"))
					message = calc(dtMolecules,dtResults,add + "ESP",canon, receivedConfiguration);

				if(message != null && (message.contains("The ESP method does not work with") || message.contains("HomoLumo")) )
					message = calc(dtMolecules,dtResults,add + AM1,canon,receivedConfiguration);

				if(message == null) break;

				if(!message.contains("This part of MOPAC is fragile")) {
					if(addedOptionalKeyword) continue; // we will still try at least one time
				}

				out.println("failed to calculate: " + message);

				if(dtResults.getRowsSize() != mol+1) // three attempts were used - row need to be added
					dtResults.addRow();
				dtResults.getCurrentRow().setError(message);

				break; // no chance to continue
			}

			while(dtResults.getRowsSize() > (mol+1)) // only the last one is correct one; two or more rows can be added by multiple processing and be wrong
				dtResults.deleteRow(mol);

		}

		return dtResults;
	}

	private String calc(DataTable dtMolecules, DataTable dtResults, String keyword, boolean canon, DescriptorsAbstractConfiguration conf){
		try{
			calculate(dtMolecules,dtResults,keyword,canon, conf);
		}
		catch(InterruptedException ee)
		{
			return "Failed by timeout.";
		}
		catch(Exception ee)
		{
			String message = getErrorMessage();
			String mess ="null";
			if(ee != null)mess = ee.getMessage();

			if(mess.contains("FINAL HEAT OF FORMATION")) // could not find anything!
				return message +" -- " + ERRORSTANDARD;
			else
				return (message + " -- " + mess).replace("Not found: CARTESIAN COORDINATES - ****", "").trim();
		}

		return null;
	}

	private String getErrorMessage() {

		String line ="";

		try{
			BufferedReader reader = getAliasedBufferedReader(getOutputFile());
			if(reader == null)return "";
			line = skipTo("*     Error and normal termination messages reported in this calculation     *", reader);
			read(reader,1);
			line = read(reader,1);
			line = line.replace('*', ' ').trim();
			reader.close();
		}catch(IOException e){
		}
		return line;
	}

	private void calculate(DataTable dtMolecules, DataTable dtResults, String keyword, boolean canon, DescriptorsAbstractConfiguration conf) throws Exception{

		int TIMEOUTMOPAC = conf.getTimeoutInSeconds();

		File f = new File(getExeFile());
		String[] commands = {"export MOPAC_LICENSE="+getAliasedPath()+";","export OMP_NUM_THREADS="+THREADS+";",getExeFile(), getAliasedFileName(FILENAME)};
		String[] env = {"LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-10.1/compat/:/usr/local/cuda-11.0/compat/:"+f.getParent()+"/:.","MOPAC_LICENSE="+f.getParent()+"/"};

		if(gpuCard != NO_GPU) {
			env = OCHEMUtils.append(env, "CUDA_VISIBLE_DEVICES=" + gpuCard);
			env = OCHEMUtils.append(env, "CUDA_DEVICE_ORDER=PCI_BUS_ID");
		}

		String keywords = getKeywords() + " T=" + (TIMEOUTMOPAC*9/10)+ " " + keyword;

		if(keyword.equals(AM1))
			keywords = keywords.replace(PM7, "");

		keywords = keywords.trim();

		prepareMOPIN(dtMolecules, keywords, canon);
		try {
			executeBinaryBash(commands, env, getOutputFile(), TIMEOUTMOPAC);
		}catch(UserFriendlyException e) {
		}

		parseMopacOutput(getOutputFile(),dtResults, conf);

		if(dtResults.getColumnIndex("TotalEnergy") != -1 && Double.isNaN((Double)dtResults.getValue("TotalEnergy")))
			throw new IOException("TotalEnergy is NaN");

		if(dtResults.getColumnIndex("MOPAC_TOTAL_ENERGY") != -1 && Double.isNaN((Double)dtResults.getValue("MOPAC_TOTAL_ENERGY")))
			throw new IOException("MOPAC_TOTAL_ENERGY is NaN");

	}


	String getOutputFile() throws IOException{
		return stdout;
	}

	/**
	 * Uses OpenBabel to create MOPAC input files from a specified SDF-file
	 *  using specified calculation keywords.
	 * @param canon 
	 * @throws InterruptedException 
	 * @throws Exception 
	 */
	private void prepareMOPIN(DataTable dtMolecules, String keywords, boolean canon) throws IOException, InterruptedException {

		String sdf = dtMolecules.getSDF();
		if(Various.molecule.getAtomCount(sdf) == 1) {
			String smiles = Various.molecule.convertToCanonicalSMILES(sdf);
			sdf = null;
			String atoms[]={"Cl","Br","I","F"};
			for(String a:atoms) 
				if(smiles.contains(a)) {
					sdf="[H]"+"["+a+"]";
					break;
				}
			if(sdf == null && smiles.contains("N"))
				sdf = "N(H)(H)(H)";
			//if(sdf != null)dtMolecules.setValue(0, sdf);
		}

		saveMolecules(dtMolecules, datainsdf, QSPRConstants.SDF, dtMolecules.currentRow, 1);

		String commands[]={ExecutableRunner.findExecutable("obabel"),"-h", "-xk",keywords, "-i" + "sdf", datainsdf,"-o"+"mopcrt","-O"+FILENAME};
		String can[]={ExecutableRunner.findExecutable("obabel"),"-h", "-xk",keywords,"--canonical","-i" + "sdf",datainsdf,"-o"+"mopcrt","-O" + FILENAME};

		executeBinary(canon ? can : commands, FILENAME, TIMEOUTBABEL);
	}

	// Parses MOPAC output file. Throws exception if parsing fails.
	protected void parseMopacOutput(String mopacFileName, DataTable dtResult, DescriptorsAbstractConfiguration conf) throws Exception
	{

		if(dtResult.getRowsSize() == 0)
			dtResult.columns.addAll(Arrays.asList("TotalEnergy","ElectronicEnergy","CoreRepulsion","FinalHeat","IonisationPotential",
					HOMO,LUMO,"DipolPointCharge","DipolHybrid","DipolSum","ElectronNuclear","ElectronElectron","ResonanceEnergy",
					"ExchangeEnergy","ElectronElectronRepulsion","ElectronNuclearAttraction","NuclearNuclearRepulsion","TotalElectrostaticInteraction",
					"PrincipalMomentsInertiaA","PrincipalMomentsInertiaB","PrincipalMomentsInertiaC"));

		dtResult.addRow();

		BufferedReader reader = getAliasedBufferedReader(mopacFileName);		

		String line; int filledLevels;

		try{

			line = skipTo("CARTESIAN COORDINATES", reader);
			read(reader,3);
			int numAtoms = 0;
			while ((line = read(reader).trim()).length() > 0)
				++numAtoms;  // Count until first empty line.

			dtResult.setValue("FinalHeat", parseValue("FINAL HEAT OF FORMATION", reader));
			//try {
				dtResult.setValue("TotalEnergy", parseValue("TOTAL ENERGY", reader));
			/*}catch(Exception e) {
				dtResult.setValue("TotalEnergy", parseValue("ETOT (EONE + ETWO)", reader));
				reader.close();
				reader = getAliasedBufferedReader(mopacFileName);
			}
			try{*/
				dtResult.setValue("ElectronicEnergy", parseValue("ELECTRONIC ENERGY", reader));
			/*}catch(Exception e) {
				reader.close();
				reader = getAliasedBufferedReader(mopacFileName);
				dtResult.setValue("ElectronicEnergy", parseValue("ELECTRON-NUCLEAR ATTRACTION", reader));

			}
			try {*/
				dtResult.setValue("CoreRepulsion", parseValue("CORE-CORE REPULSION", reader));
			/*}catch(Exception e) {
				reader.close();
				reader = getAliasedBufferedReader(mopacFileName);
				dtResult.setValue("CoreRepulsion", parseValue("NUCLEAR-NUCLEAR REPULSION", reader));
				reader.close();
				reader = getAliasedBufferedReader(mopacFileName);
			}*/

			addCOSMO(reader,dtResult);

			dtResult.setValue("IonisationPotential", parseValue("IONIZATION POTENTIAL", reader));

			// Filled levels
			line = skipTo("NO. OF FILLED LEVELS", reader);
			dtResult.setValue("NumberFilledLevels", filledLevels = Integer.parseInt(line.substring(36).trim()));

			addHomoLumo(reader,numAtoms,filledLevels,dtResult);

			addSuperDelocalizabilities(reader,dtResult);

			// Dipole contributions
			line = skipTo("DIPOLE", reader);
			dtResult.setValue("DipolPointCharge", Double.parseDouble(last(read(reader),1)));
			dtResult.setValue("DipolHybrid", Double.parseDouble(last(read(reader),1)));
			dtResult.setValue("DipolSum", Double.parseDouble(last(read(reader),1)));

			// Energy partitioning
			line = skipTo("***  SUMMARY OF ENERGY PARTITION", reader);
			read(reader,3); // Skip to first relevant entry.
			dtResult.setValue("ElectronNuclear", Double.parseDouble(last(read(reader),2)));
			dtResult.setValue("ElectronElectron", Double.parseDouble(last(read(reader),2)));

			line = skipTo("RESONANCE ENERGY", reader);
			dtResult.setValue("ResonanceEnergy", Double.parseDouble(last(line,2)));
			dtResult.setValue("ExchangeEnergy", Double.parseDouble(last(read(reader),2)));

			// Attraction and repulsion between electrons and nuclei.
			dtResult.setValue("ElectronElectronRepulsion", parseValueLine("ELECTRON-ELECTRON REPULSION", reader));
			dtResult.setValue("ElectronNuclearAttraction",  parseValueLine("ELECTRON-NUCLEAR ATTRACTION", reader));
			dtResult.setValue("NuclearNuclearRepulsion", parseValueLine("NUCLEAR-NUCLEAR REPULSION", reader));
			dtResult.setValue("TotalElectrostaticInteraction", parseValueLine("TOTAL ELECTROSTATIC INTERACTION", reader));

			// Moments of inertia
			line = skipTo("PRINCIPAL MOMENTS OF INERTIA", reader);
			line = read(reader,2); 

			dtResult.setValue("PrincipalMomentsInertiaA", Double.parseDouble(last(line, 5, "\\s+=?\\s*")));
			dtResult.setValue("PrincipalMomentsInertiaB", Double.parseDouble(last(line, 3, "\\s+=?\\s*")));
			dtResult.setValue("PrincipalMomentsInertiaC", Double.parseDouble(last(line, 1, "\\s+=?\\s*")));

			// Homo / lumo energies	

			Double homo = (Double)dtResult.getValue(HOMO);
			Double lumo = (Double)dtResult.getValue(LUMO);

			dtResult.setValue("HomoLumoGap", lumo - homo);
			dtResult.setValue("HomoLumoFraction",Math.abs((homo+lumo)/(lumo-homo)));
			dtResult.setValue("Electrophilicity",((homo * homo) + 2 * homo * lumo + (lumo * lumo))/(4 * (lumo - homo)));

			addAlphaPolar(reader,dtResult);

		}catch(Exception e){
			e.printStackTrace();
			System.out.println(e.getMessage());
			throw new IOException(e.getMessage());
		}finally{
			if(reader!=null)
				reader.close();
		}

	}

	private void addHomoLumo(BufferedReader reader, int numAtoms,int filledLevels, DataTable dtResult) throws IOException {
		String line = skipTo("EIGENVECTORS", reader);
		read(reader,2); // two empty lines.

		List<String> rows = new ArrayList<String>(numAtoms);  // the rows of the matrix (including relevant headers)
		String eigenvalues = "";
		while((line = read(reader)).startsWith("   Root No."))  // read blocks as long as there are blocks
		{
			read(reader,3); // skip next three lines.
			eigenvalues += read(reader);  // remember eigenvalues
			read(reader);  // skip next line, then read the block.

			int count = 0;
			while(!(line = read(reader)).isEmpty()) 
			{
				if(count == rows.size()) rows.add(line); else rows.set(count, rows.get(count) + line.substring(11));
				++count;
			}
			read(reader);  // skip second blank line
		}

		// Parse eigenvalues.
		String[] eigenvaluesSplit = split(eigenvalues);
		double reseigenvalues[] = new double[eigenvaluesSplit.length];
		for(int i = 0; i < reseigenvalues.length; ++i) reseigenvalues[i] = Double.parseDouble(eigenvaluesSplit[i]);

		dtResult.setValue(HOMO, reseigenvalues[filledLevels-1]);
		dtResult.setValue(LUMO, reseigenvalues[filledLevels]);
	}

	protected void addAlphaPolar(BufferedReader reader, DataTable dtResult) throws IOException{}

	protected void addSuperDelocalizabilities(BufferedReader reader, DataTable dtResult) throws IOException {}

	protected void addCOSMO(BufferedReader reader, DataTable dtResult) throws IOException {}


	// Reads from MOPAC file and does error handling.
	protected String read(BufferedReader reader, int count) throws IOException
	{
		String line = null;
		for(int i=0; i<count; i++)line = read(reader);
		return line;
	}

	// Reads from MOPAC file and does error handling.
	protected String read(BufferedReader reader) throws IOException
	{
		String line = reader.readLine();
		if(line == null) throw new IOException("Unexpected end-of-file in MOPAC output.");
		return line;
	}

	// Returns the first line starting with a given string modulo whitespace.
	protected String skipTo(String s, BufferedReader reader) throws IOException 
	{
		String previousLine="";
		try{
			String line = "";
			while(!(line.trim().startsWith(s))) {
				line = read(reader);
				if(line != null && line.trim().length() > 10)
					previousLine = " - " + line.trim();
			}
			return line;
		}catch(IOException e){
			throw new IOException("Not found: " + s + previousLine);
		}
	}

	// Splits the string by whitespace after inserting space before negative signs '-' to separate strings like "-7.5343-1.02329".
	private String[] split(String s)
	{
		return s.replaceAll("-", " -").trim().split("\\s+");
	}

	// Returns the i-th entry counting from the end after splitting s by the regular expression sep (dropping initial and trailing whitespace).
	private String last(String s,  int i,  String sep)
	{
		String[] split = s.trim().split(sep);
		return split[split.length - i];
	}

	// Convenience method. Like last(...), but with whitespace as separator. 
	protected String last(String s, int i) { return last(s, i, "\\s+"); }

	protected double parseValue(String header, BufferedReader reader) throws IOException{
		String line = skipTo(header, reader);
		if(line.contains("="))
			line = line.split("=")[1].trim().split("\\s+")[0];
		else {
			line = line.substring(line.indexOf(header) + header.length()).trim().split("\\s+")[0];
		}

		double val = Double.parseDouble(line);

		System.out.println(header + " "+ val);

		return val;
	}

	private double parseValueLine(String header, BufferedReader reader) throws IOException{
		String line = skipTo(header, reader);
		return Double.parseDouble(line.substring(line.indexOf(header)  +header.length()).trim().split("\\s+")[0]);
	}


}
