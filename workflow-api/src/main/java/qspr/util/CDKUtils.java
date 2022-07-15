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

import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.ReaderFactory;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class CDKUtils {
	
	static public IAtomContainer readOneMoleculeInAnyFormat(String anyFormatData) throws IOException{

		Reader r = null;
		ISimpleChemObjectReader reader = null;
		IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
		try {
			r = new StringReader(anyFormatData);

			ReaderFactory readerFactory = new ReaderFactory();
			reader = readerFactory.createReader(r);

			IAtomContainer mol;

			if (reader == null){
				SmilesParser smilesParser = new SmilesParser(DefaultChemObjectBuilder.getInstance());
				try{
					anyFormatData = anyFormatData.replaceAll("\\s+","").replaceAll("(^\\h*)|(\\h*$)","");
					mol = smilesParser.parseSmiles(anyFormatData);
				}catch (Exception ee){
					throw new IOException(ee.getMessage().replaceAll("(\\r|\\n)", ""));
				};
			}else
				mol = reader.read(builder.newAtomContainer());

			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol);
//			try{
//				AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
//			}catch(Error e){
//				System.out.println(e.getMessage());
//			}
			StructureDiagramGenerator gen = new StructureDiagramGenerator();
			gen.generateCoordinates(mol);

			String error = (String) mol.getProperty("Error");

			if(error != null && error.length() > 0) 
				throw new IOException(error);

			return mol;
		}catch(Exception e){
			System.out.println(e.getMessage());
			throw new IOException(e.getMessage());
		}

		finally{
			if(r!=null)r.close();
			if(reader!=null)reader.close();
		}

	}
}
