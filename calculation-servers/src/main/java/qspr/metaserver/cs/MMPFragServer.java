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
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.FileUtils;

import qspr.dao.Various;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.MMPFragConfiguration;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

/**
 * Server performs fragmentation according to Matched Molecular Pairs analysis needs
 * Forms the output datatable ready to be put the the MMP index
 */

public class MMPFragServer extends  DescriptorsAbstractExecutableServer
{

	final int maxMoleculeSize = 100;
	final int maxBreakablebondCount = 30;
	final int maxMmpFragSize = 10;

	String python = null;
	int numberFound = 0;
	static String ALL = "Al";

	Set<String> set = new HashSet<String>();

	public MMPFragServer(){
		supportedTaskType = QSPRConstants.MMPFrag;
	}

	@Override
	int getBatchSize(){
		return 1;
	}

	boolean addResulRow(DataTable dtResults, String scaffold, MMPFfr group, int maxSize){
		if(group.realAtoms() > maxSize) return false; // too large fragment to be added
		if( dtResults.getCurrentRow().size() ==0 || dtResults.getCurrentRow().getValue(0) == null){
			numberFound = 0;
			set.clear();
		}

		if(set.contains(scaffold))
			throw new UserFriendlyException("Scaffold is already added: " + scaffold);
		set.add(scaffold);
		dtResults.setValue(MMPFragConfiguration.SCAFFOLD_INCHI + numberFound, scaffold);
		dtResults.setValue(MMPFragConfiguration.FRAGMENT_INCHI + numberFound, group.inchiKey);
		dtResults.setValue(MMPFragConfiguration.FRAGMENT_SMILES+ numberFound, group.molecule);	
		numberFound++;
		return true;
	}

	@Override
	protected boolean isDataTableCompact(){
		return false;
	}

	@Override
	protected DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration configuration,
			int start, int batchSize) throws Exception {

				
		getResults().addRow().addAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM, dtMolecules.getRow(start).getAttachment(QSPRConstants.MOLECULE_ID_STEREOCHEM));

		String mols[] = Various.molecule.splitOrderedByChargeAndSize((String)dtMolecules.getValue(start,QSPRConstants.SDF_COLUMN));
		String kekule = Various.molecule.convertToKekuleSMILES((String)dtMolecules.getValue(start,QSPRConstants.SDF_COLUMN));
		System.out.println("got: " + kekule);

		if(kekule.contains(ALL))
			throw new IOException("Molecule has Al - skipping");

		boolean plus=false,minus=false;

		if(mols.length == 2)
			for(int i=0;i<mols.length;i++) {
				mols[i]=Various.molecule.convertToKekuleSMILES(mols[i]);
				if(mols[i].indexOf('-')!=-1)minus = true;
				else
					if(mols[i].indexOf('+')!=-1)plus = true;
			}

		int size = Various.molecule.getAtomCount((String)dtMolecules.getValue(start, 0))/2;

		if(mols.length == 2 && plus && minus) 
			for(int i=0; i<2; i++) {
				DataTable dtRes = new DataTable(false);
				dtRes.addRow();
				dtRes.setValue(QSPRConstants.SDF_COLUMN, mols[i]);
				calculateOne(dtRes, 0);
				readResults(size<maxMmpFragSize?size:maxMmpFragSize,stdout,mols[i==0?1:0]);
			}
		else {
			calculateOne(dtMolecules, start);
			readResults(size<maxMmpFragSize?size:maxMmpFragSize,stdout,null);
		}

		return getResults();
	}


	void calculateOne(DataTable dtMolecules, int start) throws Exception {

		saveMolecules(dtMolecules, datain, QSPRConstants.SMILESNOAROM, start, 1);

		String sdf = FileUtils.getFileAsString(getAliasedFileName(datain));

		if(Various.molecule.getAtomCount(sdf) > maxMoleculeSize)
			throw new IOException("Molecule is too large - should be smaller than " + maxMoleculeSize + " atoms but it has " + Various.molecule.getAtomCount(sdf));

		if(Various.molecule.getBreakableBoundCount(sdf) > maxBreakablebondCount)
			throw new IOException("Molecule is too large - should be smaller than " + maxBreakablebondCount + " atoms but it has " + Various.molecule.getBreakableBoundCount(sdf));

		String[] commandsLin = {python, getExeFile(), datain};
		runPython(commandsLin, stdout, CONDA.RDKIT, 10);
		python = commandsLin[0];
	}

	DataTable  readResults(int maxSize, String res, String ion) throws Exception{

		DataTable dtResults = getResults();

		BufferedReader buf = null;

		Set<MMPFfr> all = new HashSet<MMPFfr>();

		try {
			buf = getAliasedBufferedReader(stdout);

			//N.B.! MMP is defined as a group which is smaller than 50% of the scaffold (i.e., bigger group or other attached groups)

			for (MMPFrGroup fg : fragment(buf, ion)) 
			{
				boolean status = false;
				if (fg.fragments.size() == 2) {

					int n0 = 0, n1 =1;
					int num = fg.fragments.get(0).numAtoms - fg.fragments.get(1).numAtoms;
					if(num <0) {
						n0=1; n1=0;
					}

					// check size
					if(fg.fragments.get(1).numAtoms > fg.fragments.get(0).numAtoms/4)continue; // at most 25% of all atoms

					System.out.println(fg.fragments.get(n0).inchiKey + " " + fg.fragments.get(n1).molecule+ " " + fg.fragments.get(n0).molecule );

					status = addResulRow(dtResults, fg.fragments.get(n0).inchiKey, fg.fragments.get(n1), maxSize); // scaffold is larger group

				}else{
					String scaffold = "";
					int keySize = 0;
					MMPFfr group = null;

					if(fg.fragments.size() == 0) continue; // group is "unbrekable"

					for (MMPFfr f : fg.fragments){
						all.add(f);

						if (f.countAl() < fg.fragments.size() - 1)  // trying to find the main group with the largest number of attachment points
						{
							scaffold += f.inchiKey+ " "; // all detached groups will form scaffold
							keySize += f.realAtoms(); // number of atoms attached to the main group
						}
						else
							group = f; // fragment with largest number of detachments will be our MMP
					}

					if(group == null)
						System.out.println("nothing found!");

					if (group.realAtoms() >= 0.5*keySize) // number of atoms in the main group should at least two times smaller than those in the "scaffold"
						continue;

					System.out.println(scaffold +" "+ group.molecule);

					status=addResulRow(dtResults, scaffold.trim(), group, maxSize);
				}

				if(!status)
					System.out.println("Skipping: " + fg);

			}


		}catch(Exception e) {
			dtResults.getCurrentRow().setError(e.getMessage());
			System.out.print(e.getMessage());
			e.printStackTrace();
		}finally {
			if(buf!=null)buf.close();
		}

		return dtResults;
	}

	// The lowest level function, fragments a molecule
	// and returns a map of (possibly non-unique) fragment groups
	Set<MMPFrGroup> fragment(BufferedReader outputReader, String ion) throws Exception{
		Set<MMPFrGroup> fragGroupMap = new HashSet<MMPFrGroup>();

		String line;
		while((line = outputReader.readLine()) != null){
			line = line.replaceAll("\\*:\\d+", MMPFragServer.ALL);
			line = line.replaceAll("\\.", ",");
			line = line.replaceAll(",,", ",");
			String frags[] = line.split(","); // one group
			if(fragGroupMap.isEmpty()) { // add the whole molecule itself
				String mols[] = Various.molecule.getHydrogenFragmentsReplacedWithAl(frags[0]);
				for(String m: mols) {
					MMPFrGroup group = new MMPFrGroup();
					group.addFr(new MMPFfr(combine(m,ion)));
					group.addFr(new MMPFfr(MMPFragServer.ALL));
					fragGroupMap.add(group);

					if(ion != null) {
						group = new MMPFrGroup();
						group.addFr(new MMPFfr(m));
						group.addFr(new MMPFfr(ion));
						fragGroupMap.add(group);
					}
				}

			}
			MMPFrGroup group = new MMPFrGroup();

			for(int i=2; i< frags.length;i++)
				if(frags[i].length()>0) 
					group.addFr(new MMPFfr(combine(frags[i],ion)));
			fragGroupMap.add(group);
		}
		return fragGroupMap;
	}

	boolean findPair(String frag, String ion) {
		if(ion == null)return false;
		if( (frag.indexOf('+')!=-1 && ion.indexOf('-')!=-1) ||
				(frag.indexOf('-')!=-1 && ion.indexOf('+')!=-1)
				)return true;
		return false;

	}

	String combine(String frag, String ion) {
		return findPair(frag, ion)?frag + "." + ion:frag;
	}
}



class MMPFfr implements Comparable<MMPFfr>{
	public String inchiKey;
	public String molecule;
	public int numAtoms;

	public MMPFfr(String m) throws Exception {
		molecule = m.equals(MMPFragServer.ALL)?"[Al][H]":Various.molecule.convertToCanonicalSMILES(m);
		inchiKey = Various.molecule.getInChIKeyNoStero(molecule);
		numAtoms = Various.molecule.getAtomCount(molecule);

		//System.out.println(molecule+ " " + inchiKey + " " + numAtoms);

	}

	public int realAtoms() {
		return numAtoms - countAl();
	}

	public int countAl() {
		return molecule.length() - molecule.replaceAll(MMPFragServer.ALL, "A").length();
	}

	@Override
	public int hashCode() {
		return inchiKey.hashCode();
	}

	@Override
	public boolean equals(Object obj) {
		return compareTo((MMPFfr)obj) == 0;
	}

	@Override
	public int compareTo(MMPFfr f){
		if (f.numAtoms - numAtoms != 0)
			return f.numAtoms- numAtoms;
		else
			return inchiKey.compareTo(f.inchiKey);
	}
}

class MMPFrGroup{

	public List<MMPFfr> fragments = new ArrayList<MMPFfr>();

	public MMPFrGroup(){
	}

	public void addFr(MMPFfr f) {
		fragments.add(f);
		Collections.sort(fragments);
	}

	@Override
	public int hashCode(){
		return toString().hashCode();
	}

	@Override
	public String toString() {
		String s="";
		for(MMPFfr fr: fragments) {
			s += " " + fr.molecule; 
		}
		return s;
	}
}




