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

package qspr.depiction;

import java.io.BufferedOutputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import org.openscience.cdk.depict.Depiction;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import qspr.util.CDKUtils;

public class MoleculeDepictionCDK extends MoleculeDepiction {

	@Override
	public byte[] getDefaultImp() throws IOException {
		String format = toStringConfig();
		
		ByteArrayOutputStream baos = null;
		BufferedOutputStream bos = null;
		try {

			baos = new ByteArrayOutputStream();
			bos = new BufferedOutputStream(baos);
//			System.out.println(format);
			IAtomContainer mol = CDKUtils.readOneMoleculeInAnyFormat(moleculeData);
			if (format.contains("H_off")) {
				AtomContainerManipulator.suppressHydrogens(mol);
			} else {
				AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol);
			}

			double w = 150, h = 150;
			try{

				int p1 = format.indexOf("w");
				int p2 = format.indexOf(",h");
				int p3 = format.indexOf(",",p2+2);

				w= Double.parseDouble(format.substring(p1+1,p2));
				h= Double.parseDouble(format.substring(p2+2,p3));
			}catch(Exception e){
//				TODO: probably do something?
			}

			DepictionGenerator generator = new DepictionGenerator().withSize(w, h).withFillToFit()
					.withAtomColors();

			Depiction depiction = generator.depict(mol);

			depiction.writeTo(Depiction.PNG_FMT,bos);

			byte image[] = baos.toByteArray();
			return image;

		} catch (Exception e) {
			e.printStackTrace();
		}finally{
			try{
				baos.close();
				bos.close();
			}catch(Exception e){}
		}

		return null;
	}

}
