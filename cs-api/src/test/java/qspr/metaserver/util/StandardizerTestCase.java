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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.exceptions.UserFriendlyException;

import junit.framework.TestCase;
import qspr.OCHEMConfiguration;
import qspr.dao.ChemInfEngine;
import qspr.dao.Various;
import qspr.metaserver.util.MoleculeStandartizer.MolFormat;
import qspr.workflow.utils.QSPRConstants;

@RunWith(Parameterized.class)
public class StandardizerTestCase extends TestCase {

	private MoleculeStandartizer standardizer;
	private ChemInfEngine currentStandardizerOption;
	private String toDesalt;
	private String toDearomatize;
	private String randomSMILES [] = {
			"FB6(F)N1C(C)=C(C)C(C)=C1C=C5C(C)=C(C)C(CC4=C(C)C(C)=C3C=C2N(=C(C)C(C)=C2C)B(F)(F)N34)=N56",
			"CC1(C)CC2=[N]3C1=CC1=CC=C4C=C5C=CC6=[N]5[Zn]3(N14)N1C(C=CC1=C6)=C2",
			"[H]C1=C([H])C2=C([H])C3=[N@@]4([H])C(=C([H])C5=C([H])C([H])=C6C([H])=C7C([H])=C([H])C8=[N@@]7([H])[Zn@]4(N56)N2C1=C8[H])C(C([H])([H])[H])(C([H])([H])[H])C3([H])[H]"

	};

	private String[] toNormalize = {
			"Cc1ccccc1[N+](=O)[O-]",
			"Cc1ccccc1N(=O)=O",
	};
	private String[] toNeutralize = {
			"O=C(C)Oc1cc(CC[N-2])ccc1C(=O)[O-1]", // with anions
			"O=C(C)Oc1cc(CC[N+](C)(C)(C))ccc1C(=O)O", // with quaternary nitrogen
	};

	public static String canonize(String smiles) throws CDKException {
		SmilesParser   sp  = new SmilesParser(SilentChemObjectBuilder.getInstance());
		IAtomContainer mol = sp.parseSmiles(smiles);
		SmilesGenerator smigen = new SmilesGenerator(SmiFlavor.Isomeric | SmiFlavor.Canonical);
		return smigen.create(mol);
	}

	public StandardizerTestCase(ChemInfEngine engine) throws Exception {
		super();

		this.currentStandardizerOption = engine;
		Various.molecule = Various.getCheminfImpl(engine);

		this.standardizer = MoleculeStandartizer.getInstance(currentStandardizerOption,null); // libraries are assumed to be already uploaded

		this.toDesalt = Various.molecule.convertToFormat(
				"[N+]1=C2C=C(C=CC2=CC2C=CC(N)=CC1=2)N.S(=O)(=O)([O-])[O-]", 
				QSPRConstants.SDF
				);
		this.toDearomatize = Various.molecule.convertToFormat(
				"n1c2cc(N)ccc2cc3ccc(N)cc13", 
				QSPRConstants.SDF
				);

		System.out.println("Intializing standardizer test for: " + engine.toString());
	}

	@Parameterized.Parameters
	public static Collection<Object[]> input() {
		Various.getCheminfImpl(ChemInfEngine.CDK);

		ArrayList<Object[]> engines = new ArrayList<>(Arrays.asList(new Object[][] {
			new Object[] {ChemInfEngine.CDK}
		}));

		if (OCHEMConfiguration.chemaxonLicenceAvailable()) {
			try {
				Various.getCheminfImpl(ChemInfEngine.CHEMAXON);
			} catch (UserFriendlyException e) {
				System.err.println(ChemInfEngine.CHEMAXON + " is enabled, but no license was found.");
				throw new CriticalException(e);
			}
			engines.add(new Object[] {ChemInfEngine.CHEMAXON});
		}

		return engines;
	}

	@Override
	protected void setUp() throws Exception {
		standardizer = MoleculeStandartizer.getInstance(currentStandardizerOption,null); // libraries are assumed to be already uploaded
	}

	@Test
	public void testStandardize() {
		standardizer.setCleanStructure();
		standardizer.setStandardize();
		standardizer.outFormat = MolFormat.SMILES;
		try {
			String finishedProduct = standardizer.doStandartization(toNormalize[0]);
			assertEquals(canonize("Cc1ccccc1N(=O)=O"), canonize(finishedProduct));

			standardizer.setDeConvertNO2();
			finishedProduct = standardizer.doStandartization(toNormalize[0]);
			assertEquals(canonize("Cc1ccccc1[N+](=O)[O-]"), canonize(finishedProduct.replaceAll("#8", "O")));

			finishedProduct = standardizer.doStandartization(toNormalize[1]);
			assertEquals(canonize("Cc1ccccc1[N+](=O)[O-]"), canonize(finishedProduct));
		} catch (Exception e) {
			e.printStackTrace();
			fail("Standardizer failed with exception:" + e.getMessage());
		}
	}

	@Test
	public void testDeAromatize() {
		standardizer.setCleanStructure();
		standardizer.setDearomatize();
		standardizer.outFormat = MolFormat.SMILES;
		try {
			String finishedProduct = standardizer.doStandartization(toDearomatize);
			assertEquals(canonize("N1=C2C=C(N)C=CC2=CC3=CC=C(N)C=C13"), canonize(finishedProduct));
		} catch (Exception e) {
			e.printStackTrace();
			fail("Standardizer failed with exception:" + e.getMessage());
		}
	}

	@Test
	public void testDeSalt() {
		standardizer.setCleanStructure();
		standardizer.setDesalt();
		standardizer.outFormat = MolFormat.SMILES;
		try {
			String finishedProduct = standardizer.doStandartization(toDesalt);
			// FIXME: the code below fails and needs to be checked why
			//			assertEquals(canonize("[N+]1=C2C=C(N)C=CC2=CC=3C=CC(N)=CC13"), canonize(finishedProduct));
		} catch (Exception e) {
			e.printStackTrace();
			fail("Standardizer failed with exception:" + e.getMessage());
		}
	}

	@Test
	public void testNeutralize() {
		standardizer.setCleanStructure();
		standardizer.setNeutralize();
		standardizer.outFormat = MolFormat.SMILES;
		try {
			String finishedProduct = standardizer.doStandartization(toNeutralize[0]);
			assertEquals(canonize("O=C(C)OC1=CC(CCN)=CC=C1C(=O)O"), canonize(finishedProduct));

			finishedProduct = standardizer.doStandartization(toNeutralize[1]);
			assertEquals(canonize("CC(=O)Oc1cc(CC[N+](C)(C)C)ccc1C(O)=O"), canonize(finishedProduct));
		} catch (Exception e) {
			e.printStackTrace();
			fail("Standardizer failed with exception:" + e.getMessage());
		}
	}

	/*
	 * This was adapted from the main method of qspr.metaserver.util.ChemaxonStandartizer.
	 */
	@Test
	public void testRandom() {
		standardizer.setStandardize();
		standardizer.setCleanStructure();

		List<String> sdfs = null;
		try {
			sdfs = standardizer.doStandartization(Arrays.asList(randomSMILES));
		} catch (Exception e) {
			e.printStackTrace();
			fail("Standardizer failed with exception:" + e.getMessage());
		}

		for(int i = 0; i < sdfs.size(); i++) {
			try {
				String canonSMILES = Various.molecule.convertToCanonicalSMILES(sdfs.get(i));
				//				FIXME: this needs to be the same as with CHEMAXON, but maybe not possible
				//				assertEquals(randomSMILES_expected[i], canonSMILES.trim());
				System.out.println(Various.molecule.convertToCanonicalSMILES(randomSMILES[i]) +"\n" + Various.molecule.getFormula(sdfs.get(i)) + "\t" + 
						Various.molecule.getInChiKey(sdfs.get(i)) + " : " + canonSMILES);
			} catch (Exception e) {
				e.printStackTrace();
				fail("Failed to get formula for molecule:\n" + sdfs.get(i));
			}
		}
	}

}
