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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.UnsupportedEncodingException;
import java.net.MalformedURLException;
import java.net.URLEncoder;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeoutException;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.entities.Molecule;
import qspr.entities.MoleculeName;
import qspr.entities.ValidatedFact;
import qspr.util.ncbi.eSearch.ESearchResult;
import qspr.util.ncbi.eSearch.Id;
import qspr.util.ncbi.eSummary.ESummaryResult;
import qspr.util.ncbi.eSummary.Item;
import cern.colt.Timer;

public class NCBI_Utility
{
	private static transient final Logger logger = LogManager.getLogger(NCBI_Utility.class);
	private static final int SYNONYMLIST = 5;
	private static final int CANONICALSMILE = 11;

	private static final String ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pccompound&term=";
	private static final String ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pccompound&id=";
	private static final String STRUCTURE_URL = "https://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?disopt=DisplaySDF&cid=";

	public static Molecule getMoleculeByName(String moleculeName) throws MalformedURLException, JAXBException, UnsupportedEncodingException {

		if(!OCHEMConfiguration.allowExternalServices)return null;

		List<Molecule> fetchOneMolecule = getMoleculeByName(moleculeName, 1);
		if (fetchOneMolecule != null && fetchOneMolecule.size() > 0)
			return fetchOneMolecule.get(0);
		else
			return null;
	}

	public static List<Molecule> getMoleculeByName(String moleculeName, int numberOfResults) throws MalformedURLException, JAXBException, UnsupportedEncodingException {

		if (!OCHEMConfiguration.allowExternalServices || moleculeName == null || moleculeName.length() == 0)
			return new ArrayList<Molecule>(); // maybe return null here and handle that appropriately where it's used

		Timer t = new Timer(); t.start();
		List<Molecule> results = new ArrayList<Molecule>();
		String originalMoleculeName = moleculeName;
		try {

			String casrn = CASRN.checkCasrnSyntax(moleculeName);
			if (casrn != null)
			{
				moleculeName = casrn;
			}

			JAXBContext jContext = JAXBContext.newInstance("qspr.util.ncbi.eSearch:qspr.util.ncbi.eSummary");
			Unmarshaller unmarshaller = jContext.createUnmarshaller();

			numberOfResults = (numberOfResults < 0) ? 1 : numberOfResults;
			String url_encoded = ESEARCH_URL + URLEncoder.encode("\""+moleculeName+"\"", "UTF-8") + "[CSYN]&retmode=xml&retmax=" + numberOfResults;
			URLS idURL = new URLS(url_encoded);

			ESearchResult eSearchResult = (ESearchResult) unmarshaller.unmarshal(idURL.getURL());

			List<Id> idList = eSearchResult.getIdList().getId();
			boolean first = true;
			logger.info(url_encoded);
			// molecule not found || an error occured
			if (idList.isEmpty() ||
					(eSearchResult.getErrorList() != null && (eSearchResult.getErrorList().getPhraseNotFound().size() > 0 || eSearchResult.getErrorList().getFieldNotFound().size() > 0)))
			{
				// this means an empty molecule with this name is introduced as validated fact
				Molecule mol = new Molecule();
				MoleculeName molName = MoleculeName.get(originalMoleculeName);
				if (molName != null)
					mol.searchedBy = molName;
				saveMol2DB(mol, 0, first);
			}
			else
			{
				for (Id id : idList)
				{
					String pcID = id.getContent(); // logger.info("ID " + pcID);
					PubChemMolecule pcMol = new PubChemMolecule(pcID);
					pcMol.searchedBy = moleculeName;

					// get additional information
					URLS synURL = new URLS(ESUMMARY_URL + pcID);
					ESummaryResult eSummaryResult = (ESummaryResult) unmarshaller.unmarshal(synURL.getURL());

					List<Item> infoItem = eSummaryResult.getDocSum().getItem();
					pcMol.synonyms = getSynonymList(infoItem); // logger.info(pcMol.synonyms);
					pcMol.smiles = getSmiles(infoItem); // logger.info(pcMol.smiles);

					// get structure
					try {
						URLS strucURL = new URLS(STRUCTURE_URL + pcID);
						//NoS 25.05.12 Removed proprietary Sun PlainTextInputStream. Watch if it changes anything.
						InputStream stream = strucURL.openStream();
						String sdfString = convertStreamToString(stream);
						stream.close();
						pcMol.sdfString = sdfString; // logger.info(pcMol.sdfString);
						pcMol.inchi = getInchi(pcMol.sdfString);
					} catch (IOException e) {
						e.printStackTrace();
					} catch (InterruptedException e) {
						e.printStackTrace();
					}

					//pcMol.printInfoMessage(moleculeName);
					pcMol.searchedBy = originalMoleculeName; //Re
					Molecule mol = pcMol.convertPCMolToMolecule();
					results.add(mol);
					saveMol2DB(mol, Integer.parseInt(pcMol.cid), first);

					int synNum = Math.min(pcMol.synonyms.size(), 50);

					saveVFSynonym(pcMol, mol, synNum);
					first = false;
				}
			}

			logger.info(idList.size() + " Mols fetched");
			t.stop();
			logger.info("PubChem fetching time: "+ t.millis() +"ms.");

		}
		catch (Exception e)
		{
			// this means an empty molecule with this name is introduced as validated fact
			Molecule mol = new Molecule();
			// TODO refactor this
			//mol.id = NOTFOUNDINPUBCHEM; //Minus stopizzot

			MoleculeName molName = MoleculeName.get(originalMoleculeName);
			if (molName != null)
				mol.searchedBy = molName;
			saveMol2DB(mol, 0, true);
		}

		return results;
	}

	private static void saveVFSynonym(PubChemMolecule pcMol, Molecule mol, int synNum) {
		for (int i = 0; i < synNum; i++)
		{
			MoleculeName molName = MoleculeName.get(pcMol.synonyms.get(i));
			if (molName != null)
			{
				ValidatedFact vfsyn = ValidatedFact.getSynonym(molName, mol.mapping2);
				if (vfsyn == null)
				{
					vfsyn = new ValidatedFact();
					vfsyn.mapping = mol.mapping2;
					vfsyn.moleculename = molName;
					vfsyn.source = ValidatedFact.SOURCE_PUBCHEM;
					vfsyn.validated = ValidatedFact.SYNONYM;
					vfsyn.sourceid = Integer.parseInt(pcMol.cid);
				}
				Globals.session().saveOrUpdate(vfsyn);
			}
		}
	}

	private static void saveMol2DB(Molecule mol, int cid, boolean isFirst)
	{
		// be very careful here
		// mol.searchedBy might be null
		// before it contained an empty MoleculeName
		// maybe we have to check for null here

		ValidatedFact vf = ValidatedFact.getPrimary(mol.searchedBy);

		if (vf == null)
		{
			vf = ValidatedFact.getSynonym(mol.searchedBy, mol.mapping2);

			if (vf == null)
				vf = new ValidatedFact();

			vf.mapping = mol.mapping2;

			// be careful here mol.searchedBy might be null
			vf.moleculename = mol.searchedBy;

			vf.source = ValidatedFact.SOURCE_PUBCHEM;
			vf.sourceid = cid;
			if (isFirst)
				vf.validated = ValidatedFact.VALIDATED;

		} else
		{
			if (vf.source >= ValidatedFact.SOURCE_PUBCHEM)
			{
				if (mol.mapping2 != null)
				{
					vf.mapping = mol.mapping2;
					vf.sourceid = cid;
					vf.source = ValidatedFact.SOURCE_PUBCHEM;
				}
			}
		}

		Globals.session().saveOrUpdate(vf);
	}

	private static String[] getInchi(String sdfString) throws IOException, InterruptedException {
		return Molecule.getInChiKeys(sdfString);
	}

	private static String getSmiles(List<Item> infoItem) {
		if (null != infoItem.get(CANONICALSMILE).getContent()) {
			List<Object> canonicalSmiles = infoItem.get(CANONICALSMILE).getContent();
			if (canonicalSmiles.size() > 0)
				return canonicalSmiles.get(0).toString();
		}
		return "";
	}

	private static List<String> getSynonymList(List<Item> infoItem) {
		List<Object> content = infoItem.get(SYNONYMLIST).getContent();
		int limiter = 0;
		List<String> synonyms = new ArrayList<String>();
		for (Object items : content) {
			if (items instanceof Item && limiter++ < 50) {
				String synonym = ((Item) items).getContent().get(0).toString();
				synonyms.add(synonym);
			}
		}
		return synonyms;
	}

	public static String convertStreamToString(InputStream is) {
		/*
		 * To convert the InputStream to String we use the BufferedReader.readLine()
		 * method. We iterate until the BufferedReader return null which means
		 * there's no more data to read. Each line will appended to a StringBuilder
		 * and returned as String.
		 */
		BufferedReader reader = new BufferedReader(new InputStreamReader(is));
		StringBuilder sb = new StringBuilder();

		String line = null;
		try {
			while ((line = reader.readLine()) != null) {
				sb.append(line + "\n");
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				is.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		return sb.toString();
	}

	public static class PubChemMolecule{
		public String cid;
		public String sdfString = "";
		public String smiles = "";
		public String iupacName = "";
		public List<String> synonyms = null;
		public String[] inchi = {"",""};

		// TODO bad solution right now
		public String searchedBy = "";

		public PubChemMolecule(String id) {
			this.cid = id;
		}

		public Molecule convertPCMolToMolecule() throws IOException, TimeoutException 
		{
			Timer tt = new Timer(); tt.start();
			Molecule m = MoleculePeer.getMolecule(sdfString);
			m.searchSynonymList = synonyms;
			m.smile = smiles;
			MoleculeName molName = MoleculeName.get(searchedBy);
			if (molName != null)
				m.searchedBy = molName;
			tt.stop();
			return m;
		}

	}


	public static void main(String[] args) throws MalformedURLException, UnsupportedEncodingException, JAXBException
	{
		Globals.startAllTransactions();
		OCHEMConfiguration.allowExternalServices = true;
		Molecule m = NCBI_Utility.getMoleculeByName("atropine");
		System.out.println("everything was OK for " + m.smile);
		Globals.rollbackAllTransactions();
	}

}



