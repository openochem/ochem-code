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

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;

import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.SMILESWriter;
import org.openscience.cdk.normalize.Normalizer;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.w3c.dom.Document;

import org.openscience.cdk.graph.ConnectivityChecker;

import qspr.dao.Various;
import qspr.metaserver.CalculationServer;
import qspr.util.CDKUtils;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.SDFProcessor;

public class CDKStandartizer extends MoleculeStandartizer
{
	// Override the standardization methods here

	private IAtomContainer applyTransform(IAtomContainer mol, String file) throws Exception {
		InputStream standardization_config = CDKStandartizer.class.getClassLoader().getResourceAsStream(file);
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		Document config = db.parse(standardization_config);
		Normalizer.normalize(mol, config);
		return mol;
	}

	public IAtomContainer removeChargeWithHydrogens(IAtom atm) throws Exception {
		IAtomContainer mol = atm.getContainer();
		atm.setFormalCharge(0);
		CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(mol.getBuilder());
		IAtomType type = matcher.findMatchingAtomType(mol, atm);
		AtomTypeManipulator.configure(atm, type);
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(mol.getBuilder());
		adder.addImplicitHydrogens(mol, atm);
		return mol;
	}

	public List<String> doStandartization(List<String> sdfs)  throws Exception{
		List<String> results = new ArrayList<String>();
		for (String sdf : sdfs) {
			IAtomContainer mol = CDKUtils.readOneMoleculeInAnyFormat(sdf);

			if (cleanIt) {
				mol=cleanMolecule(mol);
			}

			mol = AtomContainerManipulator.removeHydrogens(mol);

			if (deSalt) {
				if (!ConnectivityChecker.isConnected(mol)) {
					IAtomContainerSet fragments = ConnectivityChecker.partitionIntoMolecules(mol);
					IAtomContainer largest_frag = fragments.getAtomContainer(0);
					for (IAtomContainer fragment :fragments.atomContainers()) {
						if (fragment.getAtomCount() > largest_frag.getAtomCount()) {
							largest_frag = fragment;
						}
					}
					mol = largest_frag;
				}
			}

			if (deAromatize) {
				Kekulization.kekulize(mol);
			}

			if (standardizeIt) {
				mol = applyTransform(mol, "standardize.xml");
			}

			if (neutralizeIt && AtomContainerManipulator.getTotalFormalCharge(mol) != 0) {
				for (IAtom atm : mol.atoms()) {
					int charge = atm.getFormalCharge();
					if (charge > 0) {
						if (!(atm.getSymbol().equals("N") && atm.getBondCount() == 4)) {
							// leave only quaternary nitrogens, otherwise remove the charge
							// FIXME: this should check valence rules and decide whether to change the charge or not based on that
							removeChargeWithHydrogens(atm);
						}
					} else if (charge < 0) {
						removeChargeWithHydrogens(atm);
					}
				}
			}


			if (deConvertNO2) {
				mol=applyTransform(mol, "deconvertNO.xml");
			}


			if (addExplicitHydrogens) {
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
			}

			if (cleanIt) {
				mol=cleanMolecule(mol);
			}

			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			BufferedOutputStream bos = new BufferedOutputStream(baos);
			String new_repr = null;			
			if (outFormat == MolFormat.SDF) {
				SDFWriter writer = new SDFWriter(bos);
				writer.write(mol);
				writer.close();
				new_repr = SDFProcessor.standartize(baos.toString());
			} else if (outFormat == MolFormat.SMILES) {
				SMILESWriter writer = new SMILESWriter(bos);
				writer.write(mol);
				writer.close();
				new_repr = baos.toString().trim();
			} else {
				new_repr = Various.molecule.convertToFormat(sdf, outFormat.toString());	
			}

			results.add(new_repr);
		}

		return results;
	}

	@Override
	public DataTable doStandartizationTable( DataTable datatable, CalculationServer server){
		long time = Calendar.getInstance().getTimeInMillis();

		datatable.reset();
		while (datatable.nextRow())
			if (!datatable.getCurrentRow().isError())
			{
				String sdf = (String)datatable.getValue();
				try
				{
					datatable.setValue(doStandartization(sdf));
				} 
				catch (Exception e)
				{
					datatable.getCurrentRow().setError(e.getMessage());
				}
				long time2 = Calendar.getInstance().getTimeInMillis();
				if (time2 - time > 1000)
				{
					setStatus("Standartized " + datatable.currentRow + " molecules out of " + datatable.getRowsSize());
					time = time2;
				}
			}

		return datatable;
	}

	public IAtomContainer cleanMolecule(IAtomContainer mol) throws IOException, CDKException {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
		CDKHydrogenAdder.getInstance(mol.getBuilder()).addImplicitHydrogens(mol);
		//		Kekulization.kekulize(mol);
		String smile = new SmilesGenerator(SmiFlavor.UseAromaticSymbols).create(mol);
		String smart_orig = smile;

		// ported from ChemaxonStandartizer
		smile = smile.replaceAll("\\[NH0\\]", "N");
		smile = smile.replaceAll("\\[NH0\\+\\]","[N+]");
		smile = smile.replaceAll("\\[nH0\\]", "n");
		smile = smile.replaceAll("\\[nH\\+\\]", "[n+]");
		smile = smile.replaceAll("\\[SH0\\]", "S");
		smile = smile.replaceAll("\\[SH0\\+\\]", "[S+]");
		smile = smile.replaceAll("\\[sH0\\]", "s");
		smile = smile.replaceAll("\\[sH0\\+\\]", "[s+]");
		smile = smile.replaceAll("\\[SeH0\\]", "[Se]");
		smile = smile.replaceAll("\\[B\\]", "B");
		smile = smile.replaceAll("\\[N\\]", "N");

		if (smart_orig != smile) {
			mol = CDKUtils.readOneMoleculeInAnyFormat(smile); // create new molecule if anything changed;
		}

		return mol;
	}

	public static void main(String[] args) throws Exception {
	}

}
