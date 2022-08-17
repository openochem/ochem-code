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

package qspr.dao;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeoutException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.utils.SDFProcessor;

/**
 *  Package to implement all operations with molecules
 * @author itetko
 *
 */

abstract public class ChemDAO {

	static Logger logger = LogManager.getLogger(ChemDAO.class);

	public ChemInfEngine engine;

	public static enum Display { VALID, ERROR };

	public static enum Properties{LOGP, HB, REFRACTIVITY, CHARGE}

	public static enum ION{ANION, CATION, MAX, ALL}

	static final String RADICAL  = "M  RAD";
	static final String SMARTS   = "M  MRV";
	static final String CHARGES  = "M  CHG";
	static final String END  = "M  END";

	public static final Integer PYTHONCGR = 39;

	public ChemDAO() {
		this.engine = engine();
	}

	public abstract ChemInfEngine engine();

	/**
	 *  Assuming that file is in 2D format
	 * @param type 
	 * @param _data
	 * @param type
	 * @return
	 */

	static public void aromatizeMolecules(DataTable molecules, Aromatisation aromatize, int start, int size) 
	{		
		final String ALREADY = "ALREADY";
		// Aromatize the molecules

		for(int i = start; i<start+size && i<molecules.getRowsSize();i++)
		{
			AbstractDataRow r = molecules.getRow(i);
			if(r.getAttachment(ALREADY) == null) // to avoid duplicated aromatizations
				try{
					String sdf = Various.molecule.aromatize((String) molecules.getValue(i, 0), aromatize);
					molecules.setValue(i, 0, sdf);
					r.addAttachment(ALREADY, "yes");
				}catch(Exception e){
					r.setError("Aromatization failed");
					molecules.setValue(i, 0,  QSPRConstants.ERROR);
					r.addAttachment(ALREADY, QSPRConstants.ERROR);
				}
		}

	}

	void testItself() throws IOException {
		convertToFormat("C",QSPRConstants.INCHIKEYS); // check that engine is functional
	}

	abstract public String aromatize(String molecule, Aromatisation type) throws IOException;

	public String convertToSmilesOrSmart(String molecule) throws IOException{
		if(molecule == null || molecule.length() == 0) throw new IOException("Empty molecule has been provided.");
		return convertToSmilesOrSmart(molecule,QSPRConstants.SMILES_FORMAT);
	}

	/**
	 * Attempts to convert to specified SMILES format; if fails converts to SMARTS
	 * @param molecule
	 * @param format  -- also SMILES but not default
	 * @return
	 * @throws IOException
	 */
	abstract public String convertToSmilesOrSmart(String molecule, String format) throws IOException;

	/**
	 * Used to unblock upload of failing molecules, which hangs the code - used in Upload of molecules
	 * @param molecule
	 * @return
	 * @throws IOException
	 * @throws TimeoutException
	 */

	public String convertToSDFUnblocked(String molecule)throws IOException, TimeoutException{
		Exception ee = null;
		try {
			return Various.molecule.convertToSDFUnblockedImp(molecule);
		}catch(Exception e){
			ee = e;
			try {
				for(ChemInfEngine engine :ChemInfEngine.values())
					if(!ChemDAO.ignoreEngine(engine)){
						ChemDAO mol = Various.getCheminfImpl(engine);
						if(mol != null)
							return mol.convertToSDFUnblockedImp(molecule);
					}
			}catch (Exception eee) {
			}
		}

		throw new IOException(ee);		
	}

	static public boolean ignoreEngine(ChemInfEngine engine){
		return engine == ChemInfEngine.NONE || Various.molecule.engine == engine;
	}

	abstract public String convertToSDFUnblockedImp(String molecule)throws IOException, TimeoutException;

	public double getMass(String sdfData){

		if(sdfData.toUpperCase().contains(QSPRConstants.ERROR))return 0;

		double mas = 0;
		try {
			mas = Various.molecule.getMassImp(sdfData);
		}catch(Throwable e) {
		}
		if(mas == 0){
			for(ChemInfEngine engine :ChemInfEngine.values())
				if(!ChemDAO.ignoreEngine(engine)) {
					ChemDAO mol = Various.getCheminfImpl(engine);
					if(mol != null) return mol.getMassImp(sdfData);
				}
		}
		return mas;
	}


	// FIXME: this is inefficient for large data -> we should also have an iterator version
	abstract public List<String> readSDFMolsFromFile(String filePath) throws IOException;

	public String convertToFormatFixMetal(String molecule, String format) throws IOException {
		return convertToFormat(MetalBondParserSdf.substituteMetalBondwithSingleBond(molecule),format);
	}

	/**
	 * Simple save to one of output formats
	 * @param molecule
	 * @return
	 * @throws IOException
	 */
	abstract public String convertToFormat(String molecule, String format)throws IOException;

	abstract public String guessFormat(String molecule) throws IOException;

	abstract public double getMassImp(String sdfData);

	abstract public DataTable getAllDataInSDF(InputStream inp) throws IOException;

	protected abstract String getMaxComponent(String molecule) throws IOException;

	static public String getInChIKeyNoSteroPart(String inchi){
		return inchi.equals(QSPRConstants.ERROR)?QSPRConstants.ERROR:inchi.substring(0,14);
	}

	public String getInChIKeyNoStero(String molecule) throws IOException{
		return getInChIKeyNoSteroPart(getInChiKey(molecule));
	}

	public String getInChiKey(String molecule){
		try {
			String smile = Various.molecule.getPropety(molecule, QSPRConstants.RDF_SMILES);
			if(smile != null) molecule = smile;
			molecule = MetalBondParserSdf.eliminateMetalBond(molecule);
			String inchi = null;
			try {
				inchi = convertToFormat(molecule, QSPRConstants.INCHIKEYS);
			}catch(IOException ee) {
				for(ChemInfEngine engine :ChemInfEngine.values())
					if(!ChemDAO.ignoreEngine(engine)) {
						ChemDAO mol = Various.getCheminfImpl(engine);
						if(mol != null)
							inchi = mol.convertToFormat(molecule, QSPRConstants.INCHIKEYS);
					}
			}
			return inchi == null? QSPRConstants.ERROR : inchi.substring(9);
		}catch(IOException e) {
			return QSPRConstants.ERROR;
		}
	}

	public String getInChiKeyDeSault(String molecule){
		try {
			String smile = Various.molecule.getPropety(molecule, QSPRConstants.RDF_SMILES);
			if(smile != null) molecule = smile;
			molecule = MetalBondParserSdf.eliminateMetalBond(molecule);
			molecule = getMaxComponent(molecule); // sdf
			return getInChiKey(molecule);
		}catch(IOException e) {
			return QSPRConstants.ERROR;
		}
	}

	public String getInChiKeyNoStereoDeSault(String molecule){
		return getInChIKeyNoSteroPart(getInChiKeyDeSault(molecule));
	}

	public List<String> getInChiKeys(List<String> sdfs)
	{
		List<String> results = new ArrayList<String>();
		for(String sdf:sdfs)
			results.add(Various.molecule.getInChiKey(SDFProcessor.standartize(sdf)));
		return results;
	}

	abstract public String getFormulaImpl(String molecule) throws IOException;

	public String getFormula(String sdfData) throws IOException{
		String formula = null;
		try {
			formula = Various.molecule.getFormulaImpl(sdfData);
		}catch(Throwable e) {
		}
		if(formula == null){
			for(ChemInfEngine engine :ChemInfEngine.values())
				if(!ChemDAO.ignoreEngine(engine)) {
					ChemDAO mol = Various.getCheminfImpl(engine);
					if(mol != null) return mol.getFormulaImpl(sdfData);
				}
		}
		return formula;
	}

	abstract public String convertToCanonicalName(String mol);

	public String convertToCanonicalSMILES(String mol) throws IOException{
		return convertToSmilesOrSmart(mol,QSPRConstants.SMILESUniqueNoHAromatic);
	}


	abstract public String convertToKekuleSMILES(String mol) throws IOException;

	abstract public String getAtomProperty(Properties property, String sdf) throws IOException;

	abstract public int  getAtomCount(String sdf) throws IOException;

	abstract public int  getBreakableBoundCount(String sdf) throws Exception;

	public String addAtomProperty(Properties property, String sdf) throws IOException{
		return sdf + "\n>  <" + property + ">" + "\n" + getAtomProperty(property, sdf) + "\n";
	}

	public String addAllHydrogens(String molecule) throws IOException {
		return convertToFormat(molecule, QSPRConstants.SDFNOAROM_WITHH);
	}

	private boolean requireModification(String s, String rules[]) {
		for(String rule: rules)
			if(s.contains(rule)) return true;
		return false;
	}

	private String convertAndCleanSDF(String molecule, String format, String rules[]) throws IOException {
		if(molecule == null || molecule.length() == 0) return molecule;
		String sdf = convertToFormat(molecule, format);
		if(!requireModification(sdf,rules))return sdf;
		String[] strParts = sdf.split("\\R");
		StringBuffer buf = new StringBuffer();
		for(String s : strParts)if(!requireModification(s, rules))
			buf.append(s+"\n");
		return convertToFormat(buf.toString(), format); // recalculate
	}

	public String addHydrogensAndRemoveRadicals(String molecule) throws IOException {
		String rules[] = {RADICAL};
		return convertAndCleanSDF(molecule, QSPRConstants.SDFNOAROM_WITHH, rules);
	}

	public String addHydrogensAndRemoveRadicalsAndSMARTS(String molecule) throws IOException {
		String rules[] = {RADICAL,SMARTS};
		return convertAndCleanSDF(molecule, QSPRConstants.SDFNOAROM_WITHH, rules);
	}

	public String addHydrogensAndRemoveAll(String molecule) throws IOException {
		if(molecule == null || molecule.length() == 0) return molecule;
		for(int i = 0; i < 2; i++) {
			String sdf = SDFProcessor.standartize(convertToFormat(molecule, QSPRConstants.SDFNOAROM_WITHH));
			String[] strParts = sdf.split("\\R");
			StringBuffer buf = new StringBuffer();
			for(String s : strParts)
				if(!s.startsWith("M") || s.startsWith(END))
					buf.append(s+"\n");
			molecule = buf.toString();
		}
		return molecule;
	}

	public String compareMolecules(String original, String optimized) {

		if(original == null || original.length() == 0) return "original molecules is empty";
		if(optimized == null || optimized.length() == 0) return "optimized molecules is empty";

		String mol1 = null, mol2 = null, smiles1 = null, smiles2=null, form1 = null, form2= null;

		try {
			mol1 = convertToFormat(original, QSPRConstants.INCHIKEYS).substring(9, 23); // ignoring stereo-chemistry changes
		}catch(Exception e) {}
		try {
			mol2 = convertToFormat(optimized, QSPRConstants.INCHIKEYS).substring(9, 23);
		}catch(Exception e) {}

		if(mol1 != null && mol2 != null && mol1.equals(mol2))return null; // identical

		try {
			form1 = getFormula(original);
			smiles1 = convertToSmilesOrSmart(original);
		}catch(Exception e) {}

		try {
			form2=getFormula(optimized);
			smiles2 = convertToSmilesOrSmart(optimized);
		}catch(Exception e) {}

		return ( form1 + " --> " + form2 + " sent: "  + smiles1 +  " got: " + smiles2);

	}

	abstract public String[] getHydrogenFragmentsReplacedWithAl(String molecule) throws Exception;

	public String standartiseSMIRKS(String smirks) throws IOException{
		smirks = smirks.replace(",",".");
		String pieces[] = smirks.split(">");
		smirks="";
		String smile="";
		try {
			for(int i=0;i<pieces.length;i++) {
				if(pieces[i].length() != 0)
					smirks += convertToKekuleSMILES(smile = pieces[i]);
				if(i != pieces.length-1)smirks += ">";
			}
		}catch(Exception e){
			String s = smile;
			try {
				if(!smile.contains("."))throw new IOException(); // only one piece ...
				String sms[] = smile.split(".");
				for(String ss: sms) 
					convertToKekuleSMILES(s=ss);
			}catch(Exception ee){
			}
			throw new IOException(s);
		}
		return smirks;
	}

	public String getPropety(String reaction, String property) {
		String reactions[] = reaction.split("\\n");
		for(int i =0; i<reactions.length;i++)
			if(reactions[i].contains(property))
				return reactions[i+1];

		return null;
	}

	int cnt; // counter
	int[] pre; // pre[v] = order in which dfs examines v
	int[] low; // low[v] = lowest preorder of any vertex connected to v
	int[][] adj;

	int dfs(int uu, int vv) 
	{
		int bridges = 0;
		pre[vv] = cnt++;
		low[vv] = pre[vv];
		for (int w = 0; w < adj[vv].length; w++) 
		{
			if (adj[vv][w] == 0)
				continue;

			if (pre[w] == -1) 
			{
				bridges += dfs(vv, w);
				low[vv] = Math.min(low[vv], low[w]);
				if (low[w] == pre[w]) 
				{
					if (adj[vv][w] == 1)
						bridges++;
				}
			} else if (w != uu)
				low[vv] = Math.min(low[vv], pre[w]);
		}
		return bridges;
	}

	/**
	 * Split molecule on fragments ordered by size
	 * @param data
	 * @return
	 */
	public abstract String[] splitOrderedByChargeAndSize(String data) throws IOException;


	public abstract boolean isFullyMappedSMIRKS(String reaction);

	/*
	private static String convertToCRSExternal(String rdf) throws IOException, InterruptedException {

		// external caching :)
		String crsfile = "/tmp/"+OCHEMUtils.getCrc32(rdf)+".crs";
		File f = new File(crsfile);
		if(f.exists() && f.length() > 2)return FileUtils.getFileAsString(crsfile);

		String CGRToSMILES = "/tmp/toSmiles.py";
		String SMILESToCRS = "/tmp/toCRS.py";
		if(!(new File(CGRToSMILES)).exists()) {
			String 	python = "import sys\n"
					+ "from CGRtools.files import *\n"
					+ "from io import StringIO\n"
					+ "\n"
					+ "file = sys.argv[1]\n"
					+ "\n"
					+ "with SDFRead(file) as r:\n"
					+ "    cgr = next(r) # OOP style application\n"
					+ "\n"
					+ "reactant_part, product_part = cgr.decompose()\n"
					+ "with open(sys.argv[2], 'w') as f:\n"
					+ "    print(reactant_part.__format__(format_spec={\"m\"})+\">>\"+product_part.__format__(format_spec={\"m\"}), file=f)\n"
					+ "f.close()";

			FileUtils.saveStringToFile(python, CGRToSMILES);

			python = "import sys\n"
					+ "from rdkit.Chem.CondensedGraphRxn import rdCondensedGraphRxn as CGR\n"
					+"import rdkit\n"
					+ "from rdkit import Chem\n"
					+ "\n"
					+"def __canoniseReaction (mix):\n"
					+ "    return '>'.join([__canonizeSmile(sm) for sm in mix.split('>')])\n"
					+ "\n"
					+ "def __canonizeSmile(sm):\n"
					+ "    try:\n"
					+ "        m = Chem.MolFromSmiles(sm, sanitize = True)\n"
					+ "        Chem.Kekulize(m)\n"
					+ "        return Chem.MolToSmiles(m, kekuleSmiles=True)\n"
					+ "    except:\n"
					+ "        m = Chem.MolFromSmiles(sm, sanitize = False)\n"
					+ "        m=rdkit.Chem.rdmolops.RemoveHs(m, sanitize = False)\n"
					+ "        return Chem.MolToSmiles(m, kekuleSmiles=False)\n"
					+ "\n"
					+ "file = open(sys.argv[1],mode='r')\n"
					+ "all_of_it = file.read()\n"
					+ "file.close()\n"
					+ "with open(sys.argv[2], 'w') as f:\n"
					+ "    print(CGR.CGRwriter(__canoniseReaction(all_of_it),  doRandom = False, charges = True), file=f)\n"
					+ "f.close()";

			FileUtils.saveStringToFile(python, SMILESToCRS);
		}

		String infile = "/tmp/"+OCHEMUtils.getCrc32(rdf)+".rdf";
		String oufile = "/tmp/"+OCHEMUtils.getCrc32(rdf)+".smi";
		FileUtils.saveStringToFile(rdf, infile);
		String[] commands = {OSType.getPython(PYTHONCGR),CGRToSMILES,infile,oufile};
		ProcessUtils.runBinary(commands, oufile, 3);
		commands[0] = OSType.getPython(PYTHONCGR);
		commands[1] = SMILESToCRS;
		commands[2] = oufile;
		commands[3] = crsfile;
		return ProcessUtils.runBinary(commands, crsfile, 3);
	}

	public String convertToCRS(String rdf){
		String crs = QSPRConstants.ERROR;
		try {
			logger.info("Internal RDF to CRS conversion unavailable. Attempting external...");
			crs = convertToCRSExternal(rdf);
			logger.info("External conversion success.");
		}catch(Exception ee) {
			logger.error("Failed to convert to CRS: " + rdf);
			ee.printStackTrace();
		}
		crs = crs.replaceAll("\\\\", "");
		crs = crs.replaceAll("/","");
		//System.out.println(smil + " --> " + crs + " for " + bond);
		return crs.trim();
	}

	public String CRSToSMILES(String crs) {
		int index;
		String crs1 = crs;
		// first part
		while((index=crs1.indexOf("{")) != -1) {
			String p1= crs1.substring(0, index);
			char c = crs1.charAt(index+1);
			if(c =='!')c='.';
			String p2= crs1.substring(index+4,crs1.length());
			crs1=p1+c+p2;
		}

		try {
			crs1 = convertToSmilesOrSmart(crs1,QSPRConstants.SMILESUniqueNoHAromatic);
		}catch(Exception e) {
			System.out.println("failed crs1");
		}

		// second part
		String crs2 = crs;
		while((index=crs2.indexOf("{")) != -1) {
			String p1= crs2.substring(0, index);
			char c = crs2.charAt(index+2);
			if(c =='!')c='.';
			String p2= crs2.substring(index+4,crs2.length());
			crs2=p1+c+p2;
		}

		try {
			crs2 = convertToSmilesOrSmart(crs2,QSPRConstants.SMILESUniqueNoHAromatic);
		}catch(Exception e) {
			System.out.println("failed crs2");
		}

		crs = crs1+">>"+crs2;

		System.out.println(crs);

		return crs;
	}

	 */

	public abstract List<String> getInChIComponents(String mol, ION type);
	public abstract double getTanimoto(byte[] fp_a, byte[] fp_b) throws Exception;
	public abstract double getTanimoto(String sdf_a, byte[] fp_b) throws Exception;
	public abstract byte[] getFingerprint(String sdf_a) throws Exception;
	public abstract String addPropertyToSDF(String sdf, String propertyName, String propertyValue) throws IOException;

}
