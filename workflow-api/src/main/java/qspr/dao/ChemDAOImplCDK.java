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

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.io.Serializable;
import java.io.StringReader;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.TimeoutException;

import org.apache.commons.lang.NotImplementedException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.charges.GasteigerMarsiliPartialCharges;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.inchi.InChIGenerator;
import org.openscience.cdk.inchi.InChIGeneratorFactory;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.io.DefaultChemObjectWriter;
import org.openscience.cdk.io.FormatFactory;
import org.openscience.cdk.io.Mol2Writer;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.formats.IChemFormat;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.similarity.Tanimoto;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;
import org.openscience.cdk.graph.Cycles;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.FileUtils;
import com.eadmet.utils.NumericalValueStandardizer;
import com.eadmet.utils.OSType;

import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.utils.SDFProcessor;
import qspr.util.CDKUtils;

public class ChemDAOImplCDK extends ChemDAO {

	private static transient final Logger logger = LogManager.getLogger(ChemDAOImplCDK.class);

	@Override
	public ChemInfEngine engine() {
		return ChemInfEngine.CDK;
	}

	@Override
	public String guessFormat(String anyFormatData) {
		Reader r = null;
		try {
			r = new StringReader(anyFormatData);
			FormatFactory factory = new FormatFactory();
			IChemFormat f = factory.guessFormat(r);
			return f == null ? "smiles": f.getPreferredNameExtension();
		}catch(IOException e){

		}finally{
			try{
				r.close();
			}catch(IOException e){}
		}

		return null;
	}

	@Override
	public String convertToSDFUnblockedImp(String anyFormatData) throws IOException {
		IAtomContainer mol = CDKUtils.readOneMoleculeInAnyFormat(anyFormatData);
		return convertToAnyFormat(mol,QSPRConstants.SDF);
	}

	@Override
	public double getMassImp(String anyFormatData){
		try {
			IAtomContainer molecule = CDKUtils.readOneMoleculeInAnyFormat(anyFormatData);
			return molecule == null? 0.0: NumericalValueStandardizer.getSignificantDigitsDouble(AtomContainerManipulator.getNaturalExactMass(molecule), NumericalValueStandardizer.SIGNIFICANT_DIGITS+1);
		}catch(IOException e) {
			return 0;
		}
	}

	@Override
	public DataTable getAllDataInSDF(InputStream ins) throws IOException {

		DataTable data = new DataTable(false);
		data.addColumn(QSPRConstants.MOLECULE);

		IteratingSDFReader reader = null;

		try{
			reader = new IteratingSDFReader(ins, DefaultChemObjectBuilder.getInstance());

			while (reader.hasNext()) {
				IAtomContainer molecule = (IAtomContainer)reader.next();
				AbstractDataRow row = data.addRow();
				row.setValue(0,SDFProcessor.standartize(convertToAnyFormat(molecule,QSPRConstants.SDF)));

				Map<Object,Object> properties = molecule.getProperties();

				for(Object property: properties.keySet()){
					if(property.toString().equals("cdk:Title") || property.toString().equals("cdk:Remark"))continue;
					Object value = properties.get(property);
					data.setValue(property.toString(), (Serializable)value);
				}
			}

		} 
		finally{
			reader.close();
		}

		return data;
	}

	@Override
	public String getFormulaImpl(String anyFormatData) throws IOException {
		IAtomContainer molecule = CDKUtils.readOneMoleculeInAnyFormat(anyFormatData);
		if(molecule == null) return null;
		IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(molecule);
		return MolecularFormulaManipulator.getString(formula);
	}

	@Override
	public String convertToCanonicalName(String anyFormatData) {
		return "CDK does not generate IUPAC names";
	}

	@Override
	public String convertToFormat(String anyFormatData, String format) throws IOException {
		try {
			IAtomContainer mol = CDKUtils.readOneMoleculeInAnyFormat(anyFormatData);
			return convertToAnyFormat(mol,format);
		}catch(Throwable e){
			e.printStackTrace();
			throw new IOException(e.getMessage());
		}
	}

	@Override
	public String convertToSmilesOrSmart(String anyFormatData, String format) throws IOException {
		return convertToFormat(anyFormatData,format);
	}

	private String convertToAnyFormat(IAtomContainer mol, String format) throws IOException {

		if(mol == null) return null;

		ByteArrayOutputStream baos = null;
		BufferedOutputStream bos = null;
		DefaultChemObjectWriter writer = null;
		SmilesGenerator sg;
		try{
			baos = new ByteArrayOutputStream();
			bos = new BufferedOutputStream(baos);

			switch(format){
			case QSPRConstants.MOL2:
				AtomContainerManipulator.suppressHydrogens(mol);
				writer = new Mol2Writer(bos);
				break;

			case QSPRConstants.SDF: 
				writer = new SDFWriter(bos);
				writer.getSetting("WriteAromaticBondTypes").setSetting("true");
				//				java.util.Properties customSettings = new java.util.Properties();
				//				 customSettings.setProperty(
				//				  "WriteAromaticBondTypes", "true"
				//				 );
				//				 PropertiesListener listener =
				//				   new PropertiesListener(customSettings);
				//				 writer.addChemObjectIOListener(listener);
				break;
			case QSPRConstants.SDFNOH:
				AtomContainerManipulator.suppressHydrogens(mol);
				writer = new SDFWriter(bos);
				writer.getSetting("WriteAromaticBondTypes").setSetting("true");
				break;
			case QSPRConstants.SDFNOAROM_WITHH:
				mol = dearomatize(mol);
				Kekulization.kekulize(mol);
			case QSPRConstants.SDFH:
				addImplicitHydrogens(mol);
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
				writer = new SDFWriter(bos);
				writer.getSetting("WriteAromaticBondTypes").setSetting("true");
				break;
			case QSPRConstants.SDFAROM_BASIC_NOH:
				writer = new SDFWriter(bos);
				writer.getSetting("WriteAromaticBondTypes").setSetting("true");
				AtomContainerManipulator.suppressHydrogens(mol);
				mol = aromatize(mol, Aromatisation.CDK_BASIC);
				break;
			case QSPRConstants.SDFAROM_BASIC_WITHH:
				writer = new SDFWriter(bos);
				writer.getSetting("WriteAromaticBondTypes").setSetting("true");
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
				mol = aromatize(mol, Aromatisation.CDK_BASIC);
				break;
			case QSPRConstants.SDFAROM_GENERAL_WITHH:
				writer = new SDFWriter(bos);
				writer.getSetting("WriteAromaticBondTypes").setSetting("true");
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
				mol = aromatize(mol, Aromatisation.CDK_GENERAL);
				break;
			case QSPRConstants.SMILESNOSTEREO:
				mol.setStereoElements(new ArrayList<IStereoElement>());
				addImplicitHydrogens(mol);
				sg  = new SmilesGenerator(SmiFlavor.Generic | SmiFlavor.UseAromaticSymbols);
				return sg.create(mol);
			case QSPRConstants.SMILESNOAROM:
				mol = dearomatize(mol);
				Kekulization.kekulize(mol);
				//				AtomContainerManipulator.suppressHydrogens(mol);
			case QSPRConstants.SMILES_FORMAT:
				addImplicitHydrogens(mol);
				sg  = new SmilesGenerator(SmiFlavor.Generic | SmiFlavor.UseAromaticSymbols);
				return sg.create(mol);
			case QSPRConstants.SMILESUniqueNoHAromatic:  // most importantly -- unique!
				mol = aromatize(mol, Aromatisation.CDK_GENERAL);
				addImplicitHydrogens(mol);
				AtomContainerManipulator.suppressHydrogens(mol);
				sg  = new SmilesGenerator(SmiFlavor.Generic | SmiFlavor.UseAromaticSymbols | SmiFlavor.Unique);
				return sg.create(mol);
			case QSPRConstants.SMILESH:  // can be also non-unique, but no hydrogens
				addImplicitHydrogens(mol);
				AtomContainerManipulator.suppressHydrogens(mol);
				sg  = new SmilesGenerator(SmiFlavor.Generic | SmiFlavor.UseAromaticSymbols);
				return sg.create(mol);
			case QSPRConstants.INCHIKEYS:
				//try {
				InChIGeneratorFactory factory = InChIGeneratorFactory.getInstance();
				InChIGenerator gen = factory.getInChIGenerator(mol);
				return "InChIKey=" + gen.getInchiKey();
				/*}catch(Throwable ee) { // FIX with Indigo


					Indigo indigo = new Indigo();

					StringWriter sw = new StringWriter();
					MDLV2000Writer writ = new MDLV2000Writer(sw);
					writ.write((IAtomContainer)mol);
					writ.close();
					IndigoObject m = indigo.loadMolecule(sw.toString());
					IndigoInchi indigoInchi = new IndigoInchi(indigo);
					String inchi = indigoInchi.getInchi(m);
					return "InChIKey=" + indigoInchi.getInchiKey(inchi);

				}
				 */
			default: throw new UserFriendlyException("The requested format \"" + format + "\" is not yet supported");
			}

			writer.write(mol);
			writer.close();
			String result = baos.toString();

			switch(format){
			case QSPRConstants.MOL2:
				result = result.substring(result.indexOf('\n') + 1).trim();
				break;
			case QSPRConstants.SMILESH:
			case QSPRConstants.SMILESNOAROM:
			case QSPRConstants.SMILESNOSTEREO:
			case QSPRConstants.SMILES_FORMAT:
			case QSPRConstants.SMILESUniqueNoHAromatic:
				result = result.trim();
				break;
			}

			return result;

		} catch (CDKException e) {
			throw new IOException(e.getMessage());
		}
		finally{
			try{
				bos.close();
				baos.close();
				if(writer!=null)writer.close();
			}catch(Exception e){
			}
		}		
	}

	private void addImplicitHydrogens(IAtomContainer mol) throws CDKException {
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
		adder.addImplicitHydrogens(mol);
	}

	@Override
	public int getAtomCount(String sdf) throws IOException {
		if(sdf == null || sdf.length() == 0)return 0;
		IAtomContainer mol = CDKUtils.readOneMoleculeInAnyFormat(sdf); 
		return AtomContainerManipulator.suppressHydrogens(mol).getAtomCount();
	}

	@Override
	public int getBreakableBoundCount(String sdf) throws Exception {
		IAtomContainer m = CDKUtils.readOneMoleculeInAnyFormat(sdf);

		low = new int[m.getAtomCount()];
		pre = new int[m.getAtomCount()];
		adj = new int[m.getAtomCount()][m.getAtomCount()];
		int bridges = 0; // number of bridges

		for (int v = 0; v < m.getAtomCount(); v++)
			low[v] = pre[v] = -1;

		for (int i = 0; i < m.getAtomCount() - 1; i++)
			for (int j = i + 1;  j < m.getAtomCount(); j++)
			{

				IBond mb = m.getBond(m.getAtom(i), m.getAtom(j));
				if (mb != null) adj[i][j] = adj[j][i] = mb.getOrder().numeric();
			}

		for (int v = 0; v < pre.length; v++)
			if (pre[v] == -1)
				bridges += dfs(v, v);

		return bridges;
	}

	/**
	 * Creates a new molecule instance with a new atom added to a given atom.
	 * 
	 * @param mol the molecule to add the atom to
	 * @param toAddTo the atom to attach the new atom to
	 * @param toBeAddedElement element of the added atom specified as given in the periodic table 
	 * @param bondOrder the order of the bond between the original atom and the new atom
	 * @return new molecule with the atom added
	 * @throws CDKException
	 * @throws CloneNotSupportedException
	 */
	public IAtomContainer addAtom(IAtomContainer mol, IAtom toAddTo, String toBeAddedElement, IBond.Order bondOrder) throws CDKException, CloneNotSupportedException {
		IAtomContainer nm = mol.clone();
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		IAtom a = builder.newInstance(IAtom.class, toBeAddedElement);
		nm.addAtom(a);
		IBond b = builder.newInstance(IBond.class);
		b.setOrder(bondOrder);
		b.setAtom(a, 0);
		IAtom toAddToInNew = nm.getAtom(toAddTo.getIndex());
		b.setAtom(toAddToInNew, 1);
		nm.addBond(b);
		CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(nm.getBuilder());
		IAtomType type = matcher.findMatchingAtomType(nm, a);
		AtomTypeManipulator.configure(a, type);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(nm.getBuilder());
		adder.addImplicitHydrogens(nm, a);
		adder.addImplicitHydrogens(nm, toAddToInNew);

		return nm;
	}

	@Override
	public String[] getHydrogenFragmentsReplacedWithAl(String molecule) throws IOException, CloneNotSupportedException, CDKException {
		Set<String> l = new HashSet<String>();
		IAtomContainer m = CDKUtils.readOneMoleculeInAnyFormat(molecule);
		m = AtomContainerManipulator.removeHydrogens(m); // for testing only, keep commented out - IVT: I think it is required...

		for (int i=0; i<m.getAtomCount(); i++)
		{
			IAtom ma = m.getAtom(i);
			if (ma.getImplicitHydrogenCount() == 0)
				continue;

			IAtomContainer nm = addAtom(m, ma, "Al", IBond.Order.SINGLE);

			l.add(convertToAnyFormat(nm, QSPRConstants.SMILESUniqueNoHAromatic));
		}

		System.out.println(l);
		return l.toArray(new String[l.size()]);
	}

	@Override
	public String convertToKekuleSMILES(String mol) throws IOException {
		IAtomContainer m = CDKUtils.readOneMoleculeInAnyFormat(mol);
		return convertToAnyFormat(m, QSPRConstants.SMILESNOAROM);
	}

	@Override
	public String aromatize(String molecule, Aromatisation type) throws IOException {
		IAtomContainer m = CDKUtils.readOneMoleculeInAnyFormat(molecule);
		m = aromatize(m, type);
		return convertToAnyFormat(m, QSPRConstants.SDF);
	}

	private IAtomContainer aromatize(IAtomContainer mol, Aromatisation type) {
		CycleFinder cycles = Cycles.or(Cycles.all(), Cycles.all(6));
		Aromaticity aromaticity = new Aromaticity(
				ElectronDonation.cdkAllowingExocyclic(),
				cycles);
		try {
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		} catch (CDKException e1) {
			e1.printStackTrace();
		}
		switch (type) {
		case CHEMAXON_BASIC:
		case CDK_BASIC:
			try {
				aromaticity.apply(mol);
			} catch (CDKException e) {
				logger.error(e);
				e.printStackTrace();
			}
			//			try {
			//				Aromaticity.cdkLegacy().apply(mol);
			//			} catch (CDKException e) {
			//				// TODO Auto-generated catch block
			//				e.printStackTrace();
			//			}
			break;

		case CHEMAXON_GENERAL:
		case CDK_GENERAL:
			cycles = Cycles.mcb();
			aromaticity = new Aromaticity(
					ElectronDonation.cdkAllowingExocyclic(),
					cycles);
			try {
				aromaticity.apply(mol);
			} catch (CDKException e) {
				logger.error(e);
				e.printStackTrace();
			}
			break;

		default:
			throw new UserFriendlyException("Unknown aromatization model: " + type.toString());
		}
		return mol;
	}

	private IAtomContainer dearomatize(IAtomContainer ac) {
		Iterator<IAtom> atoms = ac.atoms().iterator();
		while (atoms.hasNext())
			atoms.next().setFlag(CDKConstants.ISAROMATIC, false);
		Iterator<IBond> bonds = ac.bonds().iterator();
		while (bonds.hasNext())
			bonds.next().setFlag(CDKConstants.ISAROMATIC, false);
		return ac;
	}

	public int getFormalCharge(IAtomContainer mol) {
		int charge = 0;
		for (IAtom atm : mol.atoms()) {
			charge += atm.getFormalCharge();
		}

		return charge;
	}

	public double getMass(IAtomContainer mol) {
		WeightDescriptor mass_desc = new WeightDescriptor();
		DescriptorValue val = mass_desc.calculate(mol);
		return Double.parseDouble(val.getValue().toString());
	}

	public IAtomContainerSet getFragments(String mol) throws IOException {
		IAtomContainer m = CDKUtils.readOneMoleculeInAnyFormat(mol);
		m = AtomContainerManipulator.removeHydrogens(m);
		return ConnectivityChecker.partitionIntoMolecules(m);
	}

	@Override
	protected String getMaxComponent(String mol) throws IOException {
		IAtomContainerSet frags = getFragments(mol);
		int fragmentCount = frags.getAtomContainerCount();
		IAtomContainer largest = frags.getAtomContainer(0);
		for(int i=1; i<fragmentCount;i++) {
			IAtomContainer thisFragment = frags.getAtomContainer(i);
			double diff = largest.getAtomCount()-thisFragment.getAtomCount();
			if(diff<0) largest = thisFragment;
			if(diff != 0)continue;

			double largest_mass = getMass(largest);
			double this_mass = getMass(thisFragment);
			diff = largest_mass - this_mass;
			if(diff<0) largest = thisFragment;
			if(diff != 0)continue;

			if(getFormalCharge(largest) > getFormalCharge(thisFragment))  
				largest = thisFragment;
		}

		return convertToAnyFormat(largest, QSPRConstants.SDF);
	}

	class Fragment implements Comparable<Fragment> {
		IAtomContainer mol;
		String inchies = null;

		String getInchies() {
			if(inchies == null)try{
				inchies=convertToAnyFormat(mol, QSPRConstants.INCHIKEYS);
				inchies = inchies.substring(9, 23); // no stereo part
			}catch(Exception e) {
				inchies = ""+ mol.hashCode();
			}
			return inchies;
		}

		Fragment(IAtomContainer mol){
			this.mol = mol;
		}

		@Override
		public int compareTo(Fragment molecule) {
			int charge_a = getFormalCharge(molecule.mol);
			int charge_b = getFormalCharge(mol);
			int diff = charge_a - charge_b;  // ordering first by charge 
			if(diff == 0) diff = molecule.mol.getAtomCount() - mol.getAtomCount(); // by number of atoms
			if(diff == 0) diff = (int)((getMass(molecule.mol) - getMass(mol))*1000); // by molecular weight
			if(diff == 0) diff = molecule.getInchies().compareTo(getInchies()); // by inchies
			return diff;
		}
	}

	@Override
	public String[] splitOrderedByChargeAndSize(String data) throws IOException {
		IAtomContainerSet frags = getFragments(data);
		List<Fragment> frs = new ArrayList<Fragment>();
		for(IAtomContainer fr : frags.atomContainers())
			frs.add(new Fragment(fr));
		Collections.sort(frs);
		int i=0;
		String mols[] = new String[frags.getAtomContainerCount()];
		for(Fragment fr:frs) {
			mols[i++]= convertToAnyFormat(fr.mol, QSPRConstants.SDF);
		}
		return mols;
	}

	public int getMapIdx(IAtom atm) {
		Integer mapidx = atm.getProperty(CDKConstants.ATOM_ATOM_MAPPING);
		if (mapidx == null)
			return 0;
		return mapidx;
	}

	@Override
	public boolean isFullyMappedSMIRKS(String reaction) {
		try {
			if(!reaction.contains(">>")) {
				logger.info("no >> in reaction "  + " " + reaction);
				return false;
			}

			String molecules[] = reaction.split(">");
			IAtomContainer reagents = CDKUtils.readOneMoleculeInAnyFormat(molecules[0]);
			IAtomContainer products = CDKUtils.readOneMoleculeInAnyFormat(molecules[2]);
			Map<Integer,String> m = new HashMap<Integer,String>();

			for(IAtom a : reagents.atoms()) {
				if(a.getSymbol().equals("H"))continue;
				if(getMapIdx(a) ==0) {
					logger.info("not mapped atom in reagents: " + a.getSymbol() + " " + reaction);
					return false;
				}
				m.put(getMapIdx(a), a.getSymbol());
			}

			for(IAtom a : products.atoms()) {
				if(a.getSymbol().equals("H") || a.getSymbol().equals("U"))continue;

				if(getMapIdx(a) ==0) {
					logger.info("not mapped atom in products: " + a.getSymbol());
					return false;
				}

				// special case to report error in data
				if(!m.containsKey(getMapIdx(a))){
					logger.info("no mapped of atom in reagents: " + a.getSymbol()+ " " + getMapIdx(a) + " " + reaction);
					return false;
				}

				if(!a.getSymbol().equals(m.get(getMapIdx(a)))){
					logger.info("different symbols for the same number: " + getMapIdx(a) + " " + a.getSymbol()+ " " + m.get(getMapIdx(a))  + " " + reaction);
					return false;
				}
			}

			return true;
		} catch (Exception e) {
			e.printStackTrace();
			logger.info("exception: " + e.getMessage() + " " + reaction);
		}


		return false;
	}

	public int getHAcceptorCount(IAtom atom) {
		if (atom.getSymbol().equals("H")) {
			return -1;
		}

		IAtomContainer ac = atom.getContainer();
		// adapted from: https://github.com/cdk/cdk/blob/master/descriptor/qsarmolecular/src/main/java/org/openscience/cdk/qsar/descriptors/molecular/HBondAcceptorCountDescriptor.java
		// looking for suitable nitrogen atoms
		if (atom.getSymbol().equals("N") && atom.getFormalCharge() <= 0) {

			// excluding nitrogens that are adjacent to an oxygen
			List<IBond> bonds = ac.getConnectedBondsList(atom);
			int nPiBonds = 0;
			for (IBond bond : bonds) {
				if (bond.getOther(atom).getSymbol().equals("O")) return 0;
				if (IBond.Order.DOUBLE.equals(bond.getOrder())) nPiBonds++;
			}

			// if the nitrogen is aromatic and there are no pi bonds then it's
			// lone pair cannot accept any hydrogen bonds
			if (atom.getFlag(CDKConstants.ISAROMATIC) && nPiBonds == 0) return 0;

			return 1;
		}
		// looking for suitable oxygen atoms
		else if (atom.getSymbol().equals("O") && atom.getFormalCharge() <= 0) {
			//excluding oxygens that are adjacent to a nitrogen or to an aromatic carbon
			List<IBond> neighbours = ac.getConnectedBondsList(atom);
			for (IBond bond : neighbours) {
				IAtom neighbor = bond.getOther(atom);
				if (neighbor.getSymbol().equals("N") ||
						(neighbor.getSymbol().equals("O") &&
								neighbor.isAromatic() &&
								bond.getOrder() != IBond.Order.DOUBLE))
					return 0;
			}
			return 1;
		}

		return 0;
	}

	public int getHDonorCount(IAtom atom) {
		// checking for O and N atoms where the formal charge is >= 0
		if ((atom.getSymbol().equals("O") || atom.getSymbol().equals("N")) && atom.getFormalCharge() >= 0) {
			// implicit hydrogens
			Integer implicitH = atom.getImplicitHydrogenCount();
			if (implicitH == CDKConstants.UNSET) implicitH = 0;
			if (implicitH > 0) {
				return 1; // we skip the explicit hydrogens part cause we found implicit hydrogens
			}
			// explicit hydrogens
			IAtomContainer ac = atom.getContainer();
			List<IAtom> neighbours = ac.getConnectedAtomsList(atom);
			for (Object neighbour : neighbours) {
				if (((IAtom) neighbour).getSymbol().equals("H")) {
					return 1;
				}
			}
		}

		return 0;
	}

	@Override
	public String getAtomProperty(Properties property, String sdf) throws IOException {
		String str = "";
		try {
			IAtomContainer mol = CDKUtils.readOneMoleculeInAnyFormat(sdf);
			int count = mol.getAtomCount();

			switch(property){
			case CHARGE:
				//				based on https://www.sciencedirect.com/science/article/pii/0040402080801682?via%3Dihub
				for (int i=0; i < count; ++i) {
					GasteigerMarsiliPartialCharges charge = new GasteigerMarsiliPartialCharges();
					mol = charge.assignGasteigerMarsiliSigmaPartialCharges(mol, true);
					double increment = mol.getAtom(i).getCharge();
					str += Double.isFinite(increment) ? NumericalValueStandardizer.getSignificantDigits(increment) : "";
					if(i != count -1)str += ";";
				}break;

			case HB: 
				HBondAcceptorCountDescriptor hbacalc = new HBondAcceptorCountDescriptor();
				HBondDonorCountDescriptor hbdcalc = new HBondDonorCountDescriptor();
				for (int i=0; i < count; ++i) {
					IAtom atm = mol.getAtom(i);
					int acceptor = getHAcceptorCount(atm);
					int donor = getHDonorCount(atm);

					if (acceptor == -1 || (donor == 0 && acceptor == 0) )  str += "I";
					else
						if(acceptor * donor > 0)str += "A|D";
						else
							if(acceptor > 0)str += "A";
							else
								if(donor > 0)str += "D";
					if(i != count -1)str += ";";
				}break;

				//			case LOGP:
				//				logPPlugin logp = new logPPlugin();
				//				logp.setMolecule(mol);logp.run();
				//				for (int i=0; i < count; ++i) {
				//					double increment = logp.getAtomlogPIncrement(i);
				//					str += Double.isFinite(increment) ? NumericalValueStandardizer.getSignificantDigits(increment) : "";
				//					if(i != count -1)str += ";";
				//				}break;
				//
				//			case REFRACTIVITY:
				//				RefractivityPlugin refraction = new RefractivityPlugin();
				//				refraction.setMolecule(mol);refraction.run();
				//				for (int i=0; i < count; ++i) {
				//					double increment = refraction.getRefractivityIncrement(i);
				//					str += Double.isFinite(increment) ? NumericalValueStandardizer.getSignificantDigits(increment) : "";
				//					if(i != count -1)str += ";";
				//				}break;

			default:
				throw new UserFriendlyException("The requested property: \"" + property + "\" is not supported");
			}
		}catch(Exception ee) {
			throw new IOException(ee.getMessage());
		}

		return str;
	}

	static void  test(ChemDAO dao, String data) throws Exception, IOException, TimeoutException{

		//		String sdf = dao.convertToFormat(data, QSPRConstants.SDFH);
		//		sdf = dao.addPropertyToSDF(sdf, "Stuff", "New");
		//		System.out.println(sdf);

		//		System.out.println(dao.guessFormat(data));
		//		System.out.println(dao.getFormula(data));
		//		System.out.println(dao.getMass(data));
		//		System.out.println(dao.convertToSmilesOrSmart(data));
		//		System.out.println(dao.convertToSDFUnblocked(data));
		//
		//		InputStream stream = new ByteArrayInputStream(data.getBytes(StandardCharsets.UTF_8));		
		//		DataTable val = dao.getAllDataInSDF(stream);
		//		val.print(System.out);
		//
		// 				System.out.println(dao.convertToFormat(data, QSPRConstants.MOL2));
		// 				System.out.println(dao.convertToFormat(data, QSPRConstants.SDF));
		//                              System.out.println(dao.convertToFormat(data, QSPRConstants.SDFNOH));
		//                                 System.out.println(dao.convertToFormat(data, QSPRConstants.SDFNOAROM_WITHH));
		//                                 System.out.println(dao.convertToFormat(data, QSPRConstants.SDFH));
		//				System.out.println(dao.convertToFormat(data, QSPRConstants.SMILESNOSTEREO));
		// 				System.out.println(dao.convertToFormat(data, QSPRConstants.SMILESNOAROM));
		// 				System.out.println(dao.convertToFormat(data, QSPRConstants.SMILES_FORMAT));
		// 				System.out.println(dao.convertToFormat(data, QSPRConstants.SMILESUniqueNoHAromatic));
		//                                 System.out.println(dao.convertToFormat(data, QSPRConstants.SMILESH));
		// 				System.out.println(dao.convertToFormat(data, QSPRConstants.INCHIKEYS));
		//                 System.out.println(dao.convertToFormat(data, QSPRConstants.SDFAROM_BASIC_NOH));
		// 		System.out.println(dao.convertToFormat(data, QSPRConstants.SDFAROM_BASIC_WITHH));
		// 		System.out.println(dao.convertToFormat(data, QSPRConstants.SDFAROM_GENERAL_WITHH));
		//		
		//		System.out.println(dao.getAtomCount(data));
		//		System.out.println(dao.getBreakableBoundCount(data));
		//		dao.getHydrogenFragmentsReplacedWithAl(data);
		//		System.out.println(dao.convertToKekuleSMILES(data));
		//		System.out.println(dao.getMaxComponent("CC.CN.CC[C+]"));
		//		System.out.println(Arrays.asList(dao.splitOrderedByChargeAndSize("CC.CN.CC[C+]")));
		//		
		//		String SM = "[CH3:1][O:2][C:3]>>[CH3:1][C:2][C:3]";
		//		dao.isFullyMappedSMIRKS(SM);

		//		String RDF = 
		//				"1\n"
		//						+ "2\n"
		//						+ "3\n"
		//						+ " 26 26  0  0  0  0            999 V2000\n"
		//						+ "   -0.1467    3.1768    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "   -0.1467    2.3518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    0.5678    1.9393    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    1.2822    2.3518    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    0.5678    1.1143    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    1.2352    0.6294    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    0.9803   -0.1553    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    0.1553   -0.1553    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "   -0.3297   -0.8227    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "   -0.0997    0.6294    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    2.0198    0.8843    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    2.1913    1.6913    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    2.6329    0.3323    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    3.4175    0.5872    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    2.4614   -0.4747    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    0.6783    3.1768    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "   -0.1467    4.0018    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "   -0.9717    3.1768    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    1.4652   -0.8227    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "   -0.8843    0.8843    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    3.6725   -0.1974    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    4.2022    0.8421    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    3.1626    1.3718    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    1.6544   -0.3032    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    2.2899   -1.2817    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "    3.2684   -0.6462    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
		//						+ "  1  2  1  0  0  0  0\n"
		//						+ "  2  3  1  0  0  0  0\n"
		//						+ "  3  4  2  0  0  0  0\n"
		//						+ "  3  5  1  0  0  0  0\n"
		//						+ "  5  6  2  0  0  0  0\n"
		//						+ "  5 10  1  0  0  0  0\n"
		//						+ "  6  7  1  0  0  0  0\n"
		//						+ "  6 11  8  0  0  0  0\n"
		//						+ "  7  8  1  0  0  0  0\n"
		//						+ "  8  9  1  0  0  0  0\n"
		//						+ "  8 10  2  0  0  0  0\n"
		//						+ " 11 12  2  0  0  0  0\n"
		//						+ " 11 13  8  0  0  0  0\n"
		//						+ " 13 14  1  0  0  0  0\n"
		//						+ " 13 15  1  0  0  0  0\n"
		//						+ "  1 16  1  0  0  0  0\n"
		//						+ "  1 17  1  0  0  0  0\n"
		//						+ "  1 18  1  0  0  0  0\n"
		//						+ "  7 19  1  0  0  0  0\n"
		//						+ " 10 20  1  0  0  0  0\n"
		//						+ " 14 21  1  0  0  0  0\n"
		//						+ " 14 22  1  0  0  0  0\n"
		//						+ " 14 23  1  0  0  0  0\n"
		//						+ " 15 24  1  0  0  0  0\n"
		//						+ " 15 25  1  0  0  0  0\n"
		//						+ " 15 26  1  0  0  0  0\n"
		//						+ "M  STY  2   1 DAT   2 DAT\n"
		//						+ "M  SAL   1  2   6  11\n"
		//						+ "M  SDT   1 dynbond                                               \n"
		//						+ "M  SDD   1     2.4323    0.2169    DA    ALL  0       0  \n"
		//						+ "M  SED   1 0>1\n"
		//						+ "M  SAL   2  2  11  13\n"
		//						+ "M  SDT   2 dynbond                                               \n"
		//						+ "M  SDD   2     3.0454   -0.0802    DA    ALL  0       0  \n"
		//						+ "M  SED   2 1>0\n"
		//						+ "M  MRV SMA  11 [#6;A]\n"
		//						+ "M  END";
		//
		//		System.out.println(dao.convertToCRS(RDF));

		//		System.out.println(dao.addAtomProperty(Properties.CHARGE, data));
		//		System.out.println(dao.addAtomProperty(Properties.HB, data));
	}

	public static void main(String[] args) throws Exception {

		ChemDAOImplCDK dao = new ChemDAOImplCDK();

		// load test data
		String data = FileUtils.getFileAsString(OSType.getHome()+"temp/example.sdf");
		data = "C1=CC2=C(C=C1)OC3=CC=CC=C3C2=O";
		// 		data = "O=c1c2ccccc2oc2ccccc12";
		//		data = "C1=CC=CNC(=O)1";
		//        data = "O=C(O)[C@H]1N(C(=O)[C@H](C)C(Cc2ccccc2)S)CCC1";

		// test implementations
		System.out.println("CDK");
		test(dao, data);
		ChemDAO daoc = (ChemDAO) Class.forName("qspr.dao.ChemDAOImplChemAxon").newInstance();
		System.out.println(daoc.engine);
		test(daoc,data);
	}


	@Override
	public List<String> getInChIComponents(String mol, ION type) {
		// TODO Auto-generated method stub
		throw new NotImplementedException("cannot get inchi components in CDK, yet");
	}

	@Override
	public List<String> readSDFMolsFromFile(String filePath) throws IOException {
		File sdfFile = new File(filePath);
		IteratingSDFReader reader = new IteratingSDFReader(
				new FileInputStream(sdfFile), DefaultChemObjectBuilder.getInstance()
				);
		List<String> mols = new ArrayList<>();
		while (reader.hasNext()) {
			IAtomContainer molecule = (IAtomContainer)reader.next();
			mols.add(convertToAnyFormat(molecule, QSPRConstants.SDF));
		}
		return mols;
	}

	@Override
	public double getTanimoto(byte[] fp_a, byte[] fp_b) throws Exception {
		BitSet set_a = BitSet.valueOf(fp_a);
		BitSet set_b = BitSet.valueOf(fp_b);

		return Tanimoto.calculate(set_a, set_b);
	}

	@Override
	public double getTanimoto(String sdf_a, byte[] fp_b) throws Exception {
		return getTanimoto(getFingerprint(sdf_a), fp_b);
	}

	@Override
	public byte[] getFingerprint(String sdf_a) throws Exception {
		CircularFingerprinter fingerprinter = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP6, 2048);
		IAtomContainer mol = CDKUtils.readOneMoleculeInAnyFormat(sdf_a);
		return fingerprinter.getBitFingerprint(mol).asBitSet().toByteArray();
	}

	@Override
	public String addPropertyToSDF(String sdf, String propertyName, String propertyValue) throws IOException {
		IAtomContainer mol = CDKUtils.readOneMoleculeInAnyFormat(sdf);
		mol.setProperty(propertyName, propertyValue);
		return convertToAnyFormat(mol, QSPRConstants.SDFH);
	}

}

