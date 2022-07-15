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

package qspr.entities;

import java.io.IOException;

import javax.persistence.Entity;

import com.eadmet.utils.OCHEMUtils;

import qspr.dao.Repository;
import qspr.dao.Various;

// Class for storing original molecules
// complements and links to Molecule entity - with SDF (extSDF) molecules

@Entity
public class OriginalMolecule 
{
	public Integer id;
	
	public String data;
	
	public MoleculeFormat format;
	
	public String md5;
	
	public Molecule molecule;
	
	public OriginalMolecule(){}
	
	public OriginalMolecule(String _data) throws IOException {
		data = removeTimestamps(_data);
		md5 = OCHEMUtils.getMD5(data);
		format = Repository.molecule.getFormatByName(Various.molecule.guessFormat(data));
	}
	
	public static String removeTimestamps(String result)
	{
		// NoS 13.01.12 - commented since it messes up stereochemical smiles like "FC(F)(F)c2c(\N=C(\n1ccnc1)COCCC)ccc(Cl)c2"
		result = result.replaceFirst("JME.*", "");	// save db space
		result = result.replaceFirst("Marvin.*", ""); // save db space
		result = result.replaceFirst("OEChem.*", ""); // save db space
		result = result.replaceFirst("[^\n]*Mrv.*", ""); // save db space
		return result;
	}
	
}
