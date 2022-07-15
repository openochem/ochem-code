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

package qspr.metaserver.util;

import java.io.IOException;
import java.io.PrintWriter;
import java.security.NoSuchAlgorithmException;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.Vector;

import qspr.dao.ChemDAO.ION;
import qspr.dao.Various;
import qspr.metaserver.configurations.ValidationConfiguration.MixtureValidation;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;

/**
 * Groups records according to the molecules which are represented by unified
 * molecular IDs 
 * 
 * Used as keys in this map, i.e. size of MoleculeMappings corresponds to the number of different molecules in the dataset 
 * 
 * The unified IDs also check that molecules with the exactly same descriptors but different
 * initial IDs will be also merged List for each ID contains all rows with this
 * molecule in the dataset
 * 
 * 
 * @author itetko
 * 
 */
public class TrainingSetSplitter
{

	private PrintWriter out;

	// Random number id required to have reproducible bagging results
	private Random rand; 

	static final long serialVersionUID = 1L;

	private LabelsTable dtLabels = null;

	/**
	 * maps records to molecule.
	 * RecordID -> "MoleculeID"
	 */
	private Map<Integer, Set<Integer>> recordsToMolHashes = new HashMap<Integer, Set<Integer>>();

	/**
	 * Maps class to molecules that contain at least one record of this class
	 */
	private Map<String, Set<Integer>> classToMoleculeHash = new HashMap<String, Set<Integer>>();

	/**
	 * Maps molecules to each Class and count of instances per class/molecule
	 */
	private Map<Integer, Map<String, Integer>> moleculeHashToClassToCount = new HashMap<Integer, Map<String, Integer>>();

	/*
	 * Maps molecules to records
	 */
	private Map<Integer, TreeSet<Integer>> molHashesToRecords = new HashMap<Integer, TreeSet<Integer>>();

	private int uniqueMolecules()
	{
		return moleculeHashToClassToCount.size();
	}

	public TrainingSetSplitter(DescriptorsTable dtMolecules, LabelsTable dtLabelsOriginal, PrintWriter out, MixtureValidation mixtureValidation, int seed) throws IOException, NoSuchAlgorithmException
	{

		rand = new Random(seed); 
		
		this.out = out;
		dtLabels = dtLabelsOriginal;

		Map<String, Integer> md5s = new HashMap<String, Integer>(); // main purpose is de-duplication for molecules with identical descriptors

		out.println(" using validation by: " + mixtureValidation);

		for (int row = 0; row < dtLabels.getDataSize(); row++)
		{

			if (dtMolecules.getMoleculeUniqueId(row) == null)
				throw new IOException("Molecule ID is missing. In bagging and cross-validation, Molecule ID is required for the fold determination.");

			//This procedure performs exclusion of mixtures
			// A new mixture will get id starting from -1 in descending order
			// Otherwise indexing of a mixture is done for molids

			Set<String> mols = new HashSet<String>();

			if(mixtureValidation == null) mixtureValidation= MixtureValidation.MIXTURE; // by default we use whole structure

			switch(mixtureValidation) {
			case COMPONENT:
				MixtureAttachment ma  =  (MixtureAttachment) dtMolecules.getRawRow(row).getAttachment(QSPRConstants.MIXTURE_ATTACHMENT);
				if(ma == null) throw new UserFriendlyException("Inconsistent use of mixtures for COMPONENT validation. Should always have MixtureAttachment.");
				mols = ma.fractions.keySet(); // inchies, no SMILES
				break;
			case ALL:
				mols.addAll(Various.molecule.getInChIComponents((String)dtMolecules.getRawRow(row).getAttachment(QSPRConstants.SMILES_ATTACHMENT),ION.ALL));
			case ANIONS:
				mols.addAll(Various.molecule.getInChIComponents((String)dtMolecules.getRawRow(row).getAttachment(QSPRConstants.SMILES_ATTACHMENT),ION.ANION));
				break;
			case CATIONS:
				mols.addAll(Various.molecule.getInChIComponents((String)dtMolecules.getRawRow(row).getAttachment(QSPRConstants.SMILES_ATTACHMENT),ION.CATION));
				break;
			case MAXCOMP:
				mols.addAll(Various.molecule.getInChIComponents((String)dtMolecules.getRawRow(row).getAttachment(QSPRConstants.SMILES_ATTACHMENT),ION.MAX));
				break;
			case RECORD:
				mols.add(""+row); // just by record ID
				break;
			default:
				String md5 = dtMolecules.getMoleculeMD5(row ,false);
				//out.println("row md5: " + md5); //DEB
				Integer molid = dtMolecules.getMoleculeUniqueId(row);
				if(!md5s.containsKey(md5)) 
					md5s.put(md5, molid); // i.e., we reuse another molid, if descriptors are the the same
				mols.add(""+md5s.get(md5)); // id is not actually used, just grouping
				break;
			}

			Set<Integer> molecules=new HashSet<Integer>();

			for(String mol: mols) {
				Integer moleculeId = mol.hashCode();  // operating with hashcode since mixtures are inchies or otherwise we have integers
				TreeSet<Integer> rows = molHashesToRecords.get(moleculeId); // do we have it already?
				if (rows == null) // create new one, if we do not have it
					molHashesToRecords.put(moleculeId, rows = new TreeSet<Integer>());
				rows.add(row); // all  
				// mapping of the current row to respective molecular ids
				molecules.add(moleculeId);
			}
			recordsToMolHashes.put(row, molecules);
		}

		// Get the number of  occurrences of  each class of each property

		countOcurrencesForEachClass();

		out.println("In total records = " + dtMolecules.getDataSize() + " provided " + uniqueMolecules() + " unique compounds.");
	}


	private void countOcurrencesForEachClass() throws IOException{
		// Get the number of  occurrences of  each class of each property

		if (dtLabels.getDataSize() == 0)
			throw new UserFriendlyException("No molecules are provided: all of them"
					+ " were filtered on the previous steps, e.g. because of error in units conversion, failure in descriptor calculations, etc.");

		for (int i = 0; i < dtLabels.getDataSize(); i++)
		{
			String classLabel = dtLabels.getMoleculeLabel(i);

			for(Integer molNum :recordsToMolHashes.get(i)) { // Molecule IDs for this row; can be more than one for mixtures

				if (molNum == null)
					throw new IOException("reverseLookupMappingID does not exist for row: " + i);

				// do we have this  class in our list already?
				Set<Integer> classSet = classToMoleculeHash.get(classLabel);
				if (classSet == null)
					classSet = new HashSet<Integer>();
				classSet.add(molNum);
				classToMoleculeHash.put(classLabel, classSet);

				Map<String, Integer> classToCount = moleculeHashToClassToCount.get(molNum);
				if (classToCount == null)
					classToCount = new HashMap<String, Integer>();
				Integer count = classToCount.get(classLabel);
				classToCount.put(classLabel, count == null ? 1 : ++count);
				moleculeHashToClassToCount.put(molNum, classToCount);
			}

		}

		out.println("Total records with labels: " + dtLabels.getDataSize());
		out.println("Total unique molecules (with labels): " + molHashesToRecords.size());
		out.println("Total unique classes (with labels): " + classToMoleculeHash.size());
		out.println("Distribution of molecules per class:");
		for (Iterator<String> i = classToMoleculeHash.keySet().iterator(); i.hasNext();)
		{
			String mol = i.next();
			out.println(mol + " " + classToMoleculeHash.get(mol).size());
		}

	}

	/**
	 * Provides next class to be processed
	 * This class has minimum number of instances and has not been yet selected (number of records < bagsize)
	 * 
	 * @param classToMolecules
	 * @param processed
	 * @param bagsize
	 * @return
	 */

	private String getNextClass(Map<String, Set<Integer>> classToMolecules, Map<String, Integer> processed, int bagsize)
	{
		Iterator<String> set = classToMolecules.keySet().iterator();

		String selectedClass = null;
		int min = Integer.MAX_VALUE;

		while (set.hasNext())
		{
			String cl = set.next();
			Integer num = processed == null ? null : processed.get(cl);
			if (num != null && bagsize != 0 && num >= bagsize)
				continue; // OK, we have already selected molecules for this class
			Integer n = classToMolecules.get(cl).size();
			if (n >= min)
				continue; // This is not the class with minimum number of instances -- next time!
			min = n;
			selectedClass = cl;
		}
		return selectedClass;
	}

	private Map<Integer, Integer> createStratifiedBag(int bagSize) throws IOException
	{
		out.println("Stratified task");
		long time = Calendar.getInstance().getTimeInMillis();

		Map<Integer, Integer> trainingSetMolecules = new HashMap<Integer, Integer>();

		if (bagSize > 0 && bagSize / classToMoleculeHash.size() < 5)
			throw new IOException("Increase the dataset or number of instances per class to have on average at least 5 records per class."
					+ "\nCurrently requested bag size:" + bagSize + " and there are classes: " + classToMoleculeHash.size());

		bagSize /= classToMoleculeHash.size(); // required number of instances per class

		out.println("Preprocessing finished. in " + (Calendar.getInstance().getTimeInMillis() - time) + " ms");

		// selection is started according to class with minimal number of molecules

		String classToProcess = null;
		Set<Integer> moleculesProcessed = new HashSet<Integer>();
		// The set of molecules that we have selected to the training set
		Map<String, Integer> classSelectedCount = new HashMap<String, Integer>();
		// how many records we have selected for each class
		Map<String, Set<Integer>> classToSelectedMolecules = new HashMap<String, Set<Integer>>();

		// fill in with empty value, to avoid null pointers
		for (Iterator<String> it = classToMoleculeHash.keySet().iterator(); it.hasNext();)
		{
			classToProcess = it.next();
			classToSelectedMolecules.put(classToProcess, new HashSet<Integer>());
			classSelectedCount.put(classToProcess, Integer.valueOf(0));
		}

		while ((classToProcess = getNextClass(classToMoleculeHash, classSelectedCount, bagSize)) != null)
		{
			out.println("Analyse: " + classToProcess);

			// molecules that were selected to be inside of bag are randomly selected to fill in bagSize
			Set<Integer> moleculesForClass = classToMoleculeHash.get(classToProcess);

			if (bagSize == 0)
			{ // in case of default stratification, bagSize is set to the size of the minimum class
				bagSize = moleculesForClass.size();
				out.println("Bagsize: " + bagSize + " is set according to the size of " + classToProcess);
			}

			// Both Vectors and Sets are required for better speed
			Vector<Integer> moleculesNew = new Vector<Integer>(), moleculesSelectFrom = new Vector<Integer>();
			int moleculesSelected = 0;

			for (Iterator<Integer> mol = moleculesForClass.iterator(); mol.hasNext();)
			{
				Integer molecule = mol.next();
				if (!moleculesProcessed.contains(molecule))
					// was already selected -- do not use it
					moleculesNew.add(molecule);
				// this molecule has not been yet selected
				else
					moleculesSelected++;
			}

			// our task is to **leave at least** 33% of molecules as the bag out
			int maxSize = 2 * moleculesForClass.size() / 3;

			if (moleculesSelected >= maxSize)
				out.println("No new records will be added for " + classToProcess);

			// now we perform selection of records from moleculesNew

			out.println("Selecting " + (bagSize - classSelectedCount.get(classToProcess)) + " records using " + moleculesNew.size()
			+ " mols and having already selected " + moleculesSelected + " mols");

			while (classSelectedCount.get(classToProcess) < bagSize)
			{

				int addMolecule = moleculesNew.get((int) (Math.round(Math.floor(moleculesNew.size() * rand.nextDouble()))));
				// Random  from all with this type

				if (!moleculesProcessed.contains(addMolecule))
				{ // new molecule
					if (moleculesSelected >= maxSize)
					{
						out.println("Continue selection of " + (bagSize - classSelectedCount.get(classToProcess)) + " records using "
								+ moleculesSelectFrom.size() + " preselected molecules");
						moleculesNew = moleculesSelectFrom; // we will continue selection only from pre-selected molecules
						continue; // we will not use this and other molecules: they should be in bag*out
					}
					moleculesProcessed.add(addMolecule);
					moleculesSelectFrom.add(addMolecule);
					moleculesSelected++;
				}

				Integer mol = trainingSetMolecules.get(addMolecule);
				trainingSetMolecules.put(addMolecule, mol == null ? 1 : mol + 1);

				// now we need to increase counts of selected classes using all records of the selected molecule
				for (Iterator<Integer> record = molHashesToRecords.get(addMolecule).iterator(); record.hasNext();)
				{
					String classLabel = dtLabels.getMoleculeLabel(record.next());
					classSelectedCount.put(classLabel, classSelectedCount.get(classLabel) + 1);
					Set<Integer> v = classToSelectedMolecules.get(classLabel);
					v.add(addMolecule);
					classToSelectedMolecules.put(classLabel, v);
				}
			}

			Iterator<Map.Entry<String, Integer>> it = classSelectedCount.entrySet().iterator();
			int count = 0;
			while (it.hasNext())
			{
				Entry<String, Integer> entry = it.next();
				out.print("class: " + entry.getKey() + " instances: " + entry.getValue());
				out.print(" total molecules:" + classToMoleculeHash.get(entry.getKey()).size());
				out.println(" selected:" + classToSelectedMolecules.get(entry.getKey()).size());
				count += entry.getValue();
			}
			out.println("Total entries: " + count + " were selected in " + (Calendar.getInstance().getTimeInMillis() - time) + " ms\n");

		}

		return trainingSetMolecules;
	}

	private Map<Integer, Integer> createFoldStratified(int fold, int totalFolds) throws IOException
	{
		Map<Integer, Integer> trainingSetMolecules = new HashMap<Integer, Integer>();
		Set<Integer> eligibleMolecules = new HashSet<Integer>();

		String minimalClass = getNextClass(classToMoleculeHash, null, 0);
		out.println("Use as a stratified basis: " + minimalClass);
		int stratifiedMols = classToMoleculeHash.get(minimalClass).size(), n = 0;
		stratifiedMols = fold == totalFolds ? stratifiedMols : (stratifiedMols * (totalFolds - 1) / totalFolds);

		// we first identify molecules that will be eligible for this CV
		for (Integer mol : moleculeHashToClassToCount.keySet())
			if (fold == totalFolds || (n++ % totalFolds) != fold)
				eligibleMolecules.add(mol);

		// it may happen that no molecules from the class with minimum number of samples will be in our fold.
		// Ce la vie!
		for (String classToProcess : classToMoleculeHash.keySet())
		{
			Vector<Integer> tmp = new Vector<Integer>(); // molecules of our class
			for (Integer mol : classToMoleculeHash.get(classToProcess))
				if (eligibleMolecules.contains(mol))
					tmp.add(mol);

			for (int count = 0; count < stratifiedMols && tmp.size() > 0; count++)
			{ // we add molecules of our class by selecting them randomly
				Integer addMolecule = tmp.get((int) (Math.round(Math.floor(tmp.size() * rand.nextDouble()))));
				tmp.remove(addMolecule);
				trainingSetMolecules.put(addMolecule, 1);
			}
		}

		out.println("selected " + trainingSetMolecules.size() + " out of " + molHashesToRecords.size());

		return trainingSetMolecules;
	}

	public Map<Integer, Integer> createFold(boolean stratified, int fold, int totalFolds) throws IOException
	{
		if (stratified)
			return createFoldStratified(fold, totalFolds);
		else
			return createNormalFold(fold, totalFolds);
	}

	public Map<Integer, Integer> createBag(boolean stratified, int instances) throws IOException
	{
		if (stratified)
			return createStratifiedBag(instances);
		else
			return createClassicBag(instances);
	}

	private Map<Integer, Integer> createNormalFold(int fold, int folds)
	{
		Set<Integer> validationSetMoleculeHashes = new HashSet<Integer>();
		Map<Integer, Integer> trainingSetMolecules = new HashMap<Integer, Integer>();

		int i = 0;
		for (Integer mol : molHashesToRecords.keySet()) 
			if (fold == folds || (i++ % folds) == fold)
				validationSetMoleculeHashes.add(mol);
			else
				trainingSetMolecules.put(mol,1);

		// cleaning training set from molecules in validation sets for mixtures

		int initTraining = trainingSetMolecules.size();

		out.println("selected " + validationSetMoleculeHashes.size() + " out of " + molHashesToRecords.size()+ " for validation set and "
				+ initTraining+ " molecules in the training set"
				);

		for(Integer valid: validationSetMoleculeHashes) 
			for(Integer row: molHashesToRecords.get(valid)) // which other molecules are in these rows??
				for(Integer train:  recordsToMolHashes.get(row))
					trainingSetMolecules.remove(train);

		out.println("left in the training: " + trainingSetMolecules.size() + " out of " + initTraining + " thus removed " + (initTraining - trainingSetMolecules.size() ));

		return trainingSetMolecules;
	}

	// Bagging is done according to molecules, not to records
	// once molecule is selected n times, its all records will be added n times
	// to the training
	private Map<Integer, Integer> createClassicBag(int bagSize)
	{
		if (bagSize == 0)
			bagSize = dtLabels.getDataSize();

		Vector<Integer> moleculesSelectFrom = new Vector<Integer>();

		for (Iterator<Integer> i = molHashesToRecords.keySet().iterator(); i.hasNext();)
			moleculesSelectFrom.add(i.next());

		Map<Integer, Integer> trainingSetMolecules = new HashMap<Integer, Integer>();

		for (int count = 0; count < bagSize;)
		{
			int addMolecule = Long.valueOf(Math.round(Math.floor(molHashesToRecords.size() * rand.nextDouble()))).intValue();
			addMolecule = moleculesSelectFrom.get(addMolecule);
			Integer mol = trainingSetMolecules.get(addMolecule);
			trainingSetMolecules.put(addMolecule, mol == null ? 1 : mol + 1);
			count += molHashesToRecords.get(addMolecule).size();
		}
		return trainingSetMolecules;
	}

	public Set<Integer> getRecords(int hash)
	{
		return molHashesToRecords.get(hash);
	}

}
