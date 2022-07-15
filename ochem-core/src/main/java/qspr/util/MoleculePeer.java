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

package qspr.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeoutException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.dao.MetalBondParserSdf;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.depiction.MoleculeDepiction;
import qspr.entities.Mapping1;
import qspr.entities.Mapping2;
import qspr.entities.Molecule;
import qspr.entities.OriginalMolecule;
import qspr.entities.ValidatedFact;
import qspr.workflow.utils.SDFProcessor;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;

public class MoleculePeer 
{
	private static transient final Logger logger = LogManager.getLogger(MoleculePeer.class);

	public static Molecule getMolecule(String data) throws TimeoutException, IOException
	{
		if (data == null || data.equals(""))
			return Repository.molecule.getEmptyMolecule();

		// 1. CHECK WHETHER WE HAVE IT ALREADY IN OUR DATABASE

		OriginalMolecule original = Repository.molecule.getOriginalMoleculeByStructure(data);

		if (original != null && original.molecule != null){ // 
			if(OCHEMConfiguration.verboseMode >= 2) logger.info("Original stored molecule found in DB");
			return original.molecule;
		}

		original = new OriginalMolecule(data);

		// 2. THIS ORIGINAL IS NOT IN DB, OR WAS NOT YET CONVERTED TO SDF

		//   a) convert from arbitrary to SDF-2D
		String sdfData = Various.molecule.convertToSDFUnblocked(original.data);			

		//   c) normalize sdf, delete redundant data
		sdfData = SDFProcessor.standartize(sdfData);

		// 3. NOW CHECK IF THIS SDF already exist
		String md5 = OCHEMUtils.getMD5(sdfData);

		original.molecule = Repository.molecule.getMoleculeByMD5(md5);
		if (original.molecule != null)
		{
			Globals.session().saveOrUpdate(original);
			return original.molecule;
		}

		// 3.1. NOW CHECK IF THE MOLECULE WITH SAME DEPICTION MD5 EXISTS (12.04.10 change... we store Molecules only with unique depictions)

		// Inchi key is used to find unique depiction

		String[] inchi = Molecule.getInChiKeys(sdfData);
		String fullInchi2 = (inchi[0].length() == 32) ? inchi[0] : inchi[0] +"-"+inchi[1];

		MoleculeDepiction depiction = MoleculeDepiction.get(OCHEMConfiguration.getCheminfEngine());
		depiction.setMolecule(sdfData);
		depiction.setWidth(150);
		depiction.setHeight(150);
		if (inchi[0].length() != 14) {
			depiction.setColor("F5A9A9");
		}
		byte[] b = depiction.getImage();

		String picMd5 = OCHEMUtils.getMD5(b);

		original.molecule = Repository.molecule.getMoleculeByInvariants(picMd5, fullInchi2);

		if (original.molecule != null)
		{
			// yeah, it is here
			Globals.session().saveOrUpdate(original);
			return original.molecule;
		}

		// 3.5 MAYBE IT'S some kind of unknown earlier or EMPTY MOLECULE? CHECK IT

		double molWeight = Various.molecule.getMass(sdfData);

		Molecule molecule = null;
		if (molWeight == 0)
			molecule = Repository.molecule.getEmptyMolecule(); // If we can't calculate even weight; this is erroneous molecule
		else
		{
			// 4. FINALLY, SAVE NEW MOLECULE AND RETURN IT
			molecule = new Molecule();
			molecule.setData(sdfData);
			molecule.md5 = md5;
			molecule.molWeight = molWeight;
			molecule.updateInchi(inchi); // Not just in case, updateInchi() calls getInchi() again --> second unnecessary inchi calculation
			molecule.pictureMd5 = picMd5;
			molecule.molImage = b;
			Globals.session().saveOrUpdate(molecule);
		}

		original.molecule = molecule;

		Globals.session().saveOrUpdate(original);
		Globals.session().flush();
		return molecule;
	}

	public static Molecule getByID2(String moleculeString) throws Exception
	{
		Integer mapping2_id = Integer.valueOf(moleculeString.replaceFirst("M", ""));
		Mapping2 m = (Mapping2) Globals.session().get(Mapping2.class, mapping2_id);
		if (m != null)
		{
			return m.molecules.get(0);
		}
		else
			throw new UserFriendlyException("Error parsing a molecule: molecule with MID=\""+moleculeString+"\" not found in the database");
	}

	public static Molecule getByID1(String moleculeString) throws Exception
	{
		Long mapping1_id = Long.valueOf(moleculeString.replaceFirst("BM", ""));
		Mapping1 m = (Mapping1) Globals.session().get(Mapping1.class, mapping1_id);
		if (m != null)
		{
			return  m.molecules.get(0);
		}
		else
			throw new UserFriendlyException("Error parsing a molecule: molecule with BMID=\""+moleculeString+"\" not found in the database");

	}

	public static Molecule getByStructure(String moleculeString) throws Exception 
	{

		String[] moleculeStringPieces = moleculeString.split("-[(\\[a-zA-Z]");
		if (moleculeStringPieces.length > 50)
			throw new UserFriendlyException("Error parsing a molecule: a lot of dashes before non-numeric characters. Seems more like a name.");

		String format = null;

		try {
			format = Various.molecule.guessFormat(moleculeString);
		}catch(Exception e){
			moleculeString = MetalBondParserSdf.processBond9(moleculeString);
			format = Various.molecule.guessFormat(moleculeString);
		}

		if (format.toLowerCase().equals("name") || format.toLowerCase().startsWith("peptide"))
			throw new IOException("Error parsing a molecule: Not a structural information");

		Molecule mol = getMolecule(moleculeString);

		return mol;
	}


	@SuppressWarnings("unchecked")
	private static List<Molecule> getByValidatedFact(String moleculeString, int factType) throws Exception
	{
		List<Molecule> m = new ArrayList<Molecule>();

		Criteria c = Globals.session()
				.createCriteria(ValidatedFact.class)
				.add(Restrictions.eq("validated", factType))
				.createCriteria("moleculename")
				.add(Restrictions.eq("name", moleculeString))
				.setMaxResults(10)
				.addOrder(Order.asc("id"));
		List<ValidatedFact> results = c.list();
		if (results.size() > 0)
			for (ValidatedFact vf : results) 
			{
				if (vf.mapping == null)
					continue;
				Molecule mol = vf.mapping.molecules.get(0);
				mol.searchedBy = vf.moleculename;
				m.add(mol);
			}
		return m;
	}

	public static List<Molecule> getByValidatedSynonym(String moleculeString) throws Exception
	{
		return getByValidatedFact(moleculeString, ValidatedFact.SYNONYM);
	}

	public static List<Molecule> getByValidatedFact(String moleculeString) throws Exception
	{
		return getByValidatedFact(moleculeString, ValidatedFact.VALIDATED);
	}

	public static Molecule fetchFromString(String moleculeString, UploadContext context) throws Exception 
	{
		Molecule mol = null;

		if (context == null)
			context = new UploadContext();

		try 
		{
			mol = context.moleculeCache.get(moleculeString);

			if (mol != null)
				return mol;

			if (moleculeString.matches("M[0-9]+")) // MoleculeID
			{
				mol = getByID2(moleculeString);
				return mol;
			}

			if (moleculeString.matches("D[0-9]+")) // Depiction ID
			{
				mol = Repository.molecule.getMolecule(Long.valueOf(moleculeString.substring(1)));
				return mol;
			}

			if (moleculeString.matches("BM[0-9]+")) // MoleculeID - mapping 1
			{
				mol = getByID1(moleculeString);
				return mol;
			}

			if (moleculeString.equals("empty"))
			{
				mol = Repository.molecule.getEmptyMolecule();
				return mol;
			}

			try
			{
				mol = getByStructure(moleculeString);
				return mol;
			} catch (IOException e)
			{
				logger.error(e.getMessage());
				//Not a valid molecule structure, continue to recover it as CASRN or find it in the NCBI
			}

			//Everything failed - it must be a name... or CASRN?
			String name = CASRN.checkCasrnSyntax(moleculeString);
			if (name == null) //Not CASRN
				name = moleculeString;

			List<Molecule> resultsByVf = getByValidatedFact(moleculeString);	
			if (resultsByVf.size() > 0)
			{
				mol = resultsByVf.get(0);
				return mol;
			}

			if (!context.allowMoleculePubchemSearch)
				throw new Exception("Error parsing a molecule: Molecule not found, and search by name in PubChem is forbidden");

			if(OCHEMConfiguration.allowExternalServices){
				Globals.restartAllTransactions(true);
				mol = NCBI_Utility.getMoleculeByName(name);
				Globals.restartAllTransactions(true);
			}

			return mol;
		} catch (Exception e)
		{
			context.moleculeCache.put(moleculeString, e);
			throw e;
		} finally
		{
			if (mol != null)
				context.moleculeCache.put(moleculeString, mol);
		}	
	}

	public static Molecule getMolecule(File f) throws IOException, FileNotFoundException, IOException, TimeoutException
	{
		FileInputStream fis = new FileInputStream(f);
		int x = fis.available();
		byte b[]= new byte[x];
		fis.read(b);
		fis.close();
		return MoleculePeer.getMolecule(new String(b));
	}

	public static void main(String[] args){
		try {

			String moleculeString = "InChI=1S/C5H12O/c1-3-5(2)4-6/h5-6H,3-4H2,1-2H3";
			getByStructure(moleculeString);
			//String solvent ="Tris (tris-(hydroxymethyl) amino methane)/HCl buffer, adjusted to 5.8";
			//createSolventAttachment(solvent,new UploadContext());
		}catch (Exception e) {
		}
		System.exit(0);
	}


}
