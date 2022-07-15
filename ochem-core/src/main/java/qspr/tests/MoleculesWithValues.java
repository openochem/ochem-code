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

package qspr.tests;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeoutException;

import qspr.entities.Molecule;
import qspr.workflow.utils.QSPRConstants;
import qspr.util.MoleculePeer;

public class MoleculesWithValues
{
	List<Molecule> molecules = new ArrayList<Molecule>();
	List<Double> values = new ArrayList<Double>();

	public static MoleculesWithValues loadFromFile(InputStream is) throws NumberFormatException, IOException, TimeoutException   {

		MoleculesWithValues result = new MoleculesWithValues();
		String line = null;
		BufferedReader reader = new BufferedReader(new InputStreamReader(is));
		while ((line = reader.readLine()) != null)
		{
			if (line.contains(QSPRConstants.SMILES_FORMAT))
				continue;
			String[] data = line.replaceAll("\"", "").replaceAll("newline", "\n").split(",");
			result.molecules.add(MoleculePeer.getMolecule(data[0]));
			result.values.add(Double.parseDouble(data[1]));
		}
		reader.close();

		return result;
	}
}
