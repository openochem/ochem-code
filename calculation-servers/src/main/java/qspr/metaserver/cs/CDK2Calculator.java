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

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.charges.GasteigerMarsiliPartialCharges;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.qsar.DescriptorEngine;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.IDescriptor;
import org.openscience.cdk.qsar.IMolecularDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import qspr.metaserver.configurations.DescriptorsCDKConfiguration;
import qspr.util.CDKUtils;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.workflow.datatypes.DataTable;

public class CDK2Calculator implements AbstractCDKCalculator{

	static GasteigerMarsiliPartialCharges charges = new GasteigerMarsiliPartialCharges();
	IMolecularDescriptor descriptorCurrent = null;
	IAtomContainer moleculeCurrent = null;
	String currentDescriptorName = "";

	private List<IMolecularDescriptor> prepareDescriptors(CDKDescriptorsServer parent, DescriptorsCDKConfiguration configuration) 
	{
		List<String> classNames = DescriptorEngine.getDescriptorClassNameByPackage("org.openscience.cdk.qsar.descriptors.molecular", parent.jars);

		DescriptorEngine engine = new DescriptorEngine(classNames, DefaultChemObjectBuilder.getInstance());

		classNames = new ArrayList<String>(); // starting again to prepare from configuration

		for (String className : engine.getDescriptorClassNames()) {
			//			System.out.println(className);
			String[] dictClasses = engine.getDictionaryClass(className);
			if(dictClasses == null)continue;
			for (String dictClass : dictClasses) {
				if (configuration.descriptorTypes.contains(dictClass)
						&& !classNames.contains(className)){ // some descriptors

					if(className.contains("IPMolecularLearningDescriptor"))
						continue;				
					// belong to 2
					// classes!
					classNames.add(className);
				}
			}
		}

		List<IDescriptor> descriptors = engine
				.instantiateDescriptors(classNames);
		List<IMolecularDescriptor> desc = new ArrayList<IMolecularDescriptor>();

		for (IDescriptor d : descriptors)
			desc.add((IMolecularDescriptor) d);

		return desc;
	}

	@Override
	public DataTable calDescriptors(final CDKDescriptorsServer parent, final DataTable dtMolecules,  int start, int batchSize, DescriptorsAbstractConfiguration conf) throws IOException
	{
		final DataTable dtResults = parent.getResults();
		final List<IMolecularDescriptor> descriptors = prepareDescriptors(parent, (DescriptorsCDKConfiguration) conf);

		for(int i = start; i < start+batchSize; i++)try
		{
			dtResults.addRow();
			final int mol = i; // required to be used in run()

			Thread t = new Thread()
			{
				@Override
				public void run()
				{
					try 
					{
						String sdf = ((String)dtMolecules.getValue(mol, 0));
						//System.out.println(sdf);
						moleculeCurrent = checkAndCleanMolecule(CDKUtils.readOneMoleculeInAnyFormat(sdf));
						IAtomContainer kekule = null;

						charges.calculateCharges(moleculeCurrent);
						for (IMolecularDescriptor descriptor : descriptors)
						{
							currentDescriptorName = descriptor.getClass().toString()+": ";
							for(String n : descriptor.getDescriptorNames())
								currentDescriptorName += n+";";

							descriptorCurrent = descriptor;

							DescriptorValue value = null;

							try{
								if(currentDescriptorName.contains("IPMolecularLearningDescriptor"))throw new Exception();
								value = descriptor.calculate(moleculeCurrent);
							}catch(Exception ee){
								if(kekule == null){  // trying to work with kekulized structures instead
									kekule = moleculeCurrent.clone();
									Kekulization.kekulize(kekule); 
								}

								value = descriptor.calculate(kekule);
							}

							String[] names = value.getNames();
							IDescriptorResult result = value.getValue();

							if (result instanceof DoubleResult)
								dtResults.setValue(mol, names[0], ((DoubleResult) result).doubleValue());

							if (result instanceof IntegerResult) 
								dtResults.setValue(mol, names[0], Double.valueOf(((IntegerResult) result).intValue()));

							if (result instanceof DoubleArrayResult) 
								for (int i = 0; i < ((DoubleArrayResult) result).length(); i++)
									dtResults.setValue(mol, names[i], Double.valueOf(((DoubleArrayResult) result).get(i)));

							if (result instanceof IntegerArrayResult) 
								for (int i = 0; i < ((IntegerArrayResult) result).length(); i++)
									dtResults.setValue(mol, names[i], Double.valueOf(((IntegerArrayResult) result).get(i)));
						}

						//dtResults.printDebugOneMol(System.out);
						parent.setStatus("Processing molecule " + (mol + 1) + " out of " + dtMolecules.getRowsSize());

					}
					catch(Exception e){
						dtResults.getCurrentRow().setError(e.getMessage() + " when calculating " + currentDescriptorName);
						parent.setStatus("Molecule " + mol + " failed for calculation of " + currentDescriptorName);

						try{
							descriptorCurrent.calculate(moleculeCurrent);
						}catch(Exception ee){}

					}

				}
			};
			t.start();
			t.join(conf.getTimeoutInMinutes()*60000);
			if (t.isAlive())
			{
				t.interrupt();
				dtResults.getCurrentRow().setError("Timeout (1) of "+conf.getTimeoutInMinutes()+" minutes was reached when calculating " + currentDescriptorName);
				parent.setStatus("molecule " +mol+" failed because of timeout for calculation of "+currentDescriptorName);
			}

			// Checks for  really-really strange problem when CDK reported all zero and no error for erroneous molecule
			if (!dtResults.getCurrentRow().isError())
			{
				boolean allZero = true;

				for (int j = 0; j < dtResults.getColumnsSize(); j++)
					if (Math.abs((Double)dtResults.getValue(j)) > 1E-6)
						allZero = false;

				if (allZero)
					dtResults.getCurrentRow().setError("All calculated descriptors are zero, although no error message was reported.");
			}
		}catch(InterruptedException ee){
			dtResults.getCurrentRow().setError("Timeout (2) of "+conf.getTimeoutInMinutes()+" minutes was reached when calculating " + currentDescriptorName);
			parent.setStatus("molecule " + i +" failed because of timeout for calculation of "+currentDescriptorName);
		}

		return dtResults;
	}

	static public IAtomContainer checkAndCleanMolecule(IAtomContainer molecule)
			throws Exception {
		// add explicit H's if required
		CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
		adder.addImplicitHydrogens(molecule);
		AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);


		// do aromaticity using legacy rules
		try {
			Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.cdkAromaticSet());
			aromaticity.apply(molecule);
		} catch (CDKException e) {
			throw new CDKException("Error in aromaticity detection.");
		}

		// Do the configuration
		try {
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		} catch (CDKException e) {
			throw new CDKException("Error in atom typing: " + e.toString());
		}
		return molecule;
	}

}
