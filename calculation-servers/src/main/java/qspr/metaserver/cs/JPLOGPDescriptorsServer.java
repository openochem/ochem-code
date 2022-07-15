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

import java.util.Map;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import com.eadmet.utils.OCHEMUtils;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.cs.util.JPlogPDescriptor;
import qspr.util.CDKUtils;
import qspr.workflow.datatypes.DataTable;

public class JPLOGPDescriptorsServer extends CDKDescriptorsServer{


	public JPLOGPDescriptorsServer(){
		supportedTaskType = DescriptorsConfiguration.JPLOGP;
		repostSize = 5000;
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}

	@Override
	protected
	DataTable calculateDescriptors(final DataTable dtMolecules, DescriptorsAbstractConfiguration conf, int start,
			int batchSize) throws Exception {

		if(jars == null)
			jars = OCHEMUtils.loadJars(getExeFile());

		if(start == 0)
			setStatus("Starting processing molecules, total = " + dtMolecules.getRowsSize() + " jars = " + (jars != null));

		final DataTable dtResults = getResults();

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
						IAtomContainer molecule = CDK2Calculator.checkAndCleanMolecule(CDKUtils.readOneMoleculeInAnyFormat(sdf));

						AtomContainerManipulator.percieveAtomTypesAndConfigureUnsetProperties(molecule);
						CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(molecule.getBuilder());
						hAdder.addImplicitHydrogens(molecule);
						AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule);
						Aromaticity.cdkLegacy().apply(molecule);

						JPlogPDescriptor JPlogP = new JPlogPDescriptor();

						Map<Integer, Integer> holo = JPlogP.getMappedHologram(molecule);

						for(Integer atom: holo.keySet()) 
							dtResults.setValue("jplogp"+atom, holo.get(atom));

						setStatus("Processing molecule " + (mol + 1) + " out of " + dtMolecules.getRowsSize());

					}
					catch(Exception e){
						dtResults.getCurrentRow().setError(e.getMessage() + " when calculating molecule: " + mol);
						System.out.println(e.getMessage());
						setStatus("Molecule " + mol + " failed for calculation " + e.getMessage());
					}

				}
			};
			t.start();
			t.join(conf.getTimeoutInMinutes()*60000);
			if (t.isAlive())
			{
				t.interrupt();
				dtResults.getCurrentRow().setError("Timeout (1) of "+conf.getTimeoutInMinutes()+" minutes was reached");
				setStatus("molecule " +mol+" failed because of timeout ");
			}

		}catch(InterruptedException ee){
			dtResults.getCurrentRow().setError("Timeout (2) of "+conf.getTimeoutInMinutes()+" minutes was reached");
			setStatus("molecule " + i +" failed because of timeout");
		}

		return dtResults;
	}


}
