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

package qspr.modelling;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeoutException;

import javax.servlet.http.HttpServletRequest;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Hibernate;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.type.LongType;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.dao.Repository;
import qspr.dao.Various;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ConditionSet;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Molecule;
import qspr.entities.Tag;
import qspr.metaserver.Event;
import qspr.metaserver.util.MixtureAttachment;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;
import qspr.util.MoleculePeer;
import qspr.util.StatusTracker;
import qspr.util.UploadContext;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.parsers.SimpleParser;
import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.MAILERConstants;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;


/**
 * A class that prepares compounds based on user request - from a basket, a tag or an external file
 * This class complements the UI from select-compounds.xslt
 * 
 * @author midnighter
 */
@XmlRootElement
@SuppressWarnings("unchecked")
public class CompoundsProvider implements Serializable 
{
	private static final long serialVersionUID = 1L;
	private static transient final Logger logger = LogManager.getLogger(CompoundsProvider.class);
	public static final int NUMBER_OF_RECORDS_TO_FETCH = 25000;

	public String error;
	public String setDescription;

	public transient Basket basket = new Basket();
	private Tag tag;
	private File externalFile;

	/**
	 * Compact representation of the structures. Enables serialization of CompoundsProvider
	 */
	private ModelApplierAttachment attachment;

	public transient Event<Basket> basketLoaded = new Event<Basket>(this);
	public transient StatusTracker status = new StatusTracker(logger);

	/**
	 * When preparing a dataset, use compound references rather than SDF structures.
	 * It is useful for the clients that do not need structures to save memory and traffic.
	 */
	//	public boolean useCompoundReferences;

	/*
	 * This is useful when there are several compounds sets provided on a page and we need to distinguish between them
	 */
	public String scope = "";

	public boolean structuresLoaded = false;
	public boolean loadStructures = true;

	@XmlTransient
	public synchronized Basket getBasket()
	{
		try
		{
			if (structuresLoaded)
				return basket;

			if (attachment != null && !attachment.molecules.isEmpty())
				basket = attachment.getWorkData(loadStructures);

			if (externalFile != null)
				basket = loadFromFile(externalFile);
			else if (tag != null)
				basket = loadFromTag(tag);
			else if (basket != null && basket.getRowsSize() > 0)
			{
				if (basket.id != null)
				{
					status.set("Initializing basket entries");
					basket = Basket.getById(basket.id);
					loadBasketEntries(basket);
				} else
				{
					if (loadStructures)
						fetchStructures(basket);
				}
			}
			else
				throw new UserFriendlyException("No compounds selected");

			basketLoaded.fire(basket);
		}
		catch (Exception e)
		{
			OCHEMUtils.rethrowSafely(e);
		}
		status.set("Prepared a set of " + basket.getRowsSize() + " compounds");
		structuresLoaded = true;
		return basket;
	}

	/**
	 * Store all the structure IDs in a lightweight and easily serializable attachment object
	 */
	public void prepareForSerialization() {
		attachment = new ModelApplierAttachment();
		attachment.setWorkData(getBasket());
	}

	public boolean hasCompounds()
	{
		if (tag != null)
			return tag.getMoleculesCount() > 0;
		return (externalFile != null || basket.getRowsSize() > 0);
	}

	public long getCompoundsNum()
	{
		if (attachment != null)
			return attachment.molecules.size();
		if (basket != null && basket.getRowsSize() > 0)
			return basket.getRowsSize();
		if (tag != null)
			return tag.getMoleculesCount();
		else if (externalFile != null)
			return 1; // FIXME: Can we identify the num of compounds in a file without processing it completely?
		else
			return basket.getRowsSize();

	}

	public void setBasket(Basket basket)
	{
		structuresLoaded = false;
		this.basket = basket;
	}

	/**
	 * Stepwise loading of basket entries for large baskets
	 * @param basket
	 */
	private void loadBasketEntries(Basket basket) {
		List<Integer> epIDs = Globals.session().createCriteria(BasketEntry.class).add(Restrictions.eq("basket", basket)).setProjection(Projections.groupProperty("id")).list();
		int total = epIDs.size();

		Globals.session().evict(basket);
		basket.entries = new ArrayList<BasketEntry>();
		while (!epIDs.isEmpty())
		{
			status.set("Loaded " + basket.entries.size() + " basket entries out of " + total);
			List<Integer> batchIDs = epIDs.subList(0, Math.min(NUMBER_OF_RECORDS_TO_FETCH, epIDs.size()));
			List<BasketEntry> batchEntries = Globals.session().createCriteria(BasketEntry.class).add(Restrictions.in("id", batchIDs)).list();
			basket.entries.addAll(batchEntries);
			batchIDs.clear();
			Globals.restartAllTransactions(true);
			MemoryUtils.ensureSufficientMemory();
		}
	}

	private Basket loadFromTag(Tag tag)
	{
		Basket basket = new Basket();
		setDescription = "Tag " + tag.name;

		List<Long> ids = Globals.session().createSQLQuery(
				"select molecule_id from Molecule natural left join Mapping1 natural left join MoleculeTag where tag_id="
						+ tag.id + " group by mapping1_id order by molecule_id limit 500000").addScalar(
								"molecule_id", LongType.INSTANCE).list();
		Globals.restartAllTransactions(true);

		int i = 0;
		while (i < ids.size())
		{
			int next = i + 1000;
			if (next > ids.size())
				next = ids.size();
			List<Long> cIds = ids.subList(i, next);
			List<Molecule> mols = Globals.session().createCriteria(Molecule.class).add(
					Restrictions.in("id", cIds)).list();
			i = next;

			status.set("" + i + " compounds loaded from the tag");
			logger.info(MemoryUtils.memorySummary());
			Globals.restartAllTransactions(false);

			for (Molecule molecule : mols)
			{
				Globals.session().evict(molecule);
				molecule.molImage = null;
				ExperimentalProperty ep = new ExperimentalProperty();
				ep.molecule = molecule;
				basket.entries.add(new BasketEntry(ep));
				MemoryUtils.ensureSufficientMemory();
			}
		}

		return basket;
	}

	private Basket loadFromFile(File f) throws Exception
	{
		String name = f.getName().toLowerCase();

		if(!Globals.userSession().user.isSuperUser() && !OCHEMConfiguration.inhouseInstallation) {

			long MB = 20*1024*1024;
			long VAL = 5*MB;

			double rank = (Globals.userSession() == null || Globals.userSession().user == null)? 0: Globals.userSession().user.rank;

			rank = rank == 0? MB: VAL;

			status.set("Loading molecules from " + f.getName() + " with length " + OCHEMUtils.getSizeBytes(f.length()) +"; maximal allowed size " + OCHEMUtils.getSizeBytes(rank));

			if ( f.length() > rank) {
				String message = "Uploaded file " + f.getName() + " has length " + OCHEMUtils.getSizeBytes(f.length())+" but only file with length < " +
						OCHEMUtils.getSizeBytes(rank) + " can be uploaded. Validated users can upload " + OCHEMUtils.getSizeBytes(VAL) + " MB files. Contact " + 
						MAILERConstants.EMAIL_OCHEM + " if you need to upload larger files.";
				status.set(message);

				Mailer.postMailSafely(new Email(MAILERConstants.EMAIL_ADMIN, "User tried to upload too large file", 
						(Globals.userSession().user == null? "Anonymous " : Globals.userSession().user.login) + " " + message
						).useHTML());

				throw new UserFriendlyException(message);
			}
		}
		Basket basket = new Basket();
		if (name.contains(".xls") || name.contains(".csv") || name.contains(".smi"))
			processXLS(basket, f);
		else
			processSDF(basket, f);

		return basket;
	}

	public CompoundsProvider parseUI()
	{
		error = null;
		structuresLoaded = false;
		basket = null;
		tag = null;
		externalFile = null;
		basket = new Basket();

		try
		{
			// Check if we have really submitted the correct form
			if (ThreadScope.get().localMpRequest.getParameter("compounds-selection-form") == null)
				throw new UserFriendlyException();

			String dataSource = getParameter("data");

			if ("molecule".equals(dataSource))
			{
				ExperimentalProperty ep = new ExperimentalProperty();
				Long id = Long.valueOf(getParameter("n-molecule"));
				Molecule mol = Repository.molecule.getMolecule(id);
				ep.molecule = mol;
				basket.entries.add(new BasketEntry(ep));
				setDescription = "One molecule: M" + mol.mapping2.id;
			}
			else if ("moldata".equals(dataSource))
			{
				ExperimentalProperty ep = new ExperimentalProperty();
				String moldata = new String(getParameter("moldata"));
				Molecule mol = MoleculePeer.getMolecule(moldata);
				ep.molecule = mol;
				basket.entries.add(new BasketEntry(ep));
			}
			else if ("validationset".equals(dataSource))
			{
				basket = Basket.getBasket(Globals.userSession(), Long.valueOf(getParameter("validationsetname")));
				Hibernate.initialize(basket);
				setDescription = "Basket: " + basket.name + " (" + basket.getRowsSize() + " compounds)";
			}
			else if ("tag".equals(dataSource))
			{
				tag = (Tag) Globals.session().get(Tag.class, Long.valueOf(getParameter("tag")));
				setDescription = "Tag: " + tag.name;
			}
			else if ("externalfile".equals(dataSource))
			{
				externalFile = Globals.getUploadedFile(scope + "externalfile");
				if (externalFile == null)
					throw new UserFriendlyException("You did not select the file with compounds");
				setDescription = "Uploaded external file: " + externalFile.getName();
			}

			if (getCompoundsNum() == 0)
				throw new UserFriendlyException("No compounds selected");
		}
		catch (UserFriendlyException ue)
		{
			error = ue.getMessage();
			throw ue;
		}
		catch (Exception e)
		{
			error = e.getMessage();
			//e.printStackTrace();
			OCHEMUtils.rethrowSafely(e);
		}

		return this;
	}

	/**
	 * Fetch the structures from the database.
	 * @param basket
	 */
	private void fetchStructures(Basket basket) 
	{
		int i = 0, n = 0;
		for (BasketEntry be : basket.entries)
		{
			Molecule m = be.ep.molecule;
			if(++i >= n) {
				status.set("Loaded " + i + " structures");
				n = i*10/9;
			}
			String molData = m.getData();
			if (m.id == null)
			{
				// We got the data, but no ID. Persist the structures.
				try {
					synchronized (CompoundsProvider.class)
					{
						be.ep.molecule = MoleculePeer.getMolecule(molData);
						Globals.restartAllTransactions(true);
					}
				} catch (IOException e) {
					be.ep.molecule = Repository.molecule.getEmptyMolecule();
				}catch (TimeoutException e) {
					be.ep.molecule = Repository.molecule.getEmptyMolecule();
				}
			}
			else if (molData == null)
			{
				// We got no data. Fetch she structures from DB by ID
				be.ep.molecule = Repository.molecule.getMolecule(m.id);
			}
			else if (molData.length() < 20 && molData.matches("M[0-9]+"))
				be.ep.molecule = Repository.molecule.getMapping2(Integer.valueOf(molData.substring(1))).getMolecule();

			if (i % 500 == 0)
			{
				Globals.restartAllTransactions(true);
				MemoryUtils.ensureSufficientMemory();
			}
		}
	}

	private String getParameter(String name)
	{
		HttpServletRequest mp = ThreadScope.get().localRequest;
		return mp.getParameter(scope + name);
	}

	private void processSDF(Basket basket, File f) throws IOException
	{

		UploadContext context = new UploadContext();

		InputStream inp = null;
		DataTable allData = null;
		try{
			inp = new FileInputStream(f);
			status.set("Preloading all molecules in memory - can require a lot of time.");
			allData = Various.molecule.getAllDataInSDF(inp);
		}catch(IOException e){

		}finally{
			if(inp!=null)inp.close();
		}

		status.set("Preloading completed");

		int mixture = allData.containsColumn(QSPRConstants.MIXTURE_CONDITION) ? allData.getColumnIndex(QSPRConstants.MIXTURE_CONDITION) : -1;

		for(int i = 0 ; i < allData.getRowsSize(); i++){

			ExperimentalProperty ep = new ExperimentalProperty();

			try 
			{
				ep.molecule = MoleculePeer.getMolecule((String)allData.getValue(i, 0));
				if(allData.getRow(i).isError())ep.errorComment = allData.getRow(i).detailedStatus;
				if(mixture > 0)
					processMixture(ep, (String)allData.getValue(i, mixture), context);

				ep.other = "";
				for(int c = 1;c<allData.getColumnsSize();c++) {
					String s = ""+allData.getValue(i, c);
					s = s.replaceAll("\t","");
					ep.other += (s.length()<100?allData.getColumn(c)+":" + s + "\t\n":"");
				}

			} catch (Exception e)
			{
				System.out.println("molecule n="+i+" failed due to " + e.getMessage());
				ep.molecule = Repository.molecule.getEmptyMolecule(); // to preserve order and number of molecules
			}

			basket.entries.add(new BasketEntry(ep));

			if (i % 100 == 0)
			{
				status.set("" + i + " compounds loaded");
				logger.info("" + i + " compounds processed. Restarting transaction...");
				Globals.restartAllTransactions(true);
			}
		}
	}

	void processMixture(ExperimentalProperty ep, String mixtureMol, UploadContext context) throws Exception{
		MixtureAttachment at = ExperimentalProperty.createMixtureAttachment(mixtureMol);
		ep.molecule = MoleculePeer.fetchFromString(at.smiles(), context);
		if(at.fractions != null) {
			if(ep.conditions == null) ep.conditions = new ConditionSet();
			ep.conditions.addMixtures(at.toString());
		}
	}

	private void processXLS(Basket basket, File f) throws Exception{
		SimpleParser parser = SimpleParser.getParser(f.getName()).setSource(f);
		List<String> headerCols = parser.sheetColumns.get(0);
		int smilesIndex = -1;
		for (int i = 0; i < headerCols.size(); i++) 
			if (QSPRConstants.SMILES_FORMAT.equalsIgnoreCase(headerCols.get(i)))
			{
				smilesIndex = i;
				break;
			}
		
		if (smilesIndex == -1 && f.getName().contains(".smi"))
			smilesIndex = 0;

		if (smilesIndex == -1) 
			throw new UserFriendlyException("There is no SMILES column found in the uploaded file. Please name the column with molecules 'SMILES'");

		while (parser.hasNext())
		{
			if ((parser.currentRow + 1) % 100 == 0)
				Globals.restartAllTransactions(true);

			status.set("Processing row " + (parser.currentRow + 1) + " out of " + (parser.lastRow + 1) + " from the Excel file");

			ExperimentalProperty ep = new ExperimentalProperty();			
			Molecule mol;
			try
			{
				List<String> values = parser.next();
				String molecule = values.get(smilesIndex);

				if (molecule == null || molecule.equals(""))
					break;

				mol = MoleculePeer.getMolecule(molecule);
			} catch (Exception e) 
			{
				mol = Molecule.getStub();
				mol.error = e.getMessage();
			}

			ep.molecule = mol;
			basket.entries.add(new BasketEntry(ep));
		}
		Globals.restartAllTransactions(true);
	}

}
