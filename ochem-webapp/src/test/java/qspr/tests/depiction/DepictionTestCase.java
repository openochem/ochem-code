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

package qspr.tests.depiction;

import org.junit.Test;

import qspr.dao.ChemInfEngine;

public class DepictionTestCase extends DepictionTestCaseBase {
	
	public DepictionTestCase(ChemInfEngine engine) throws Exception {
		super(engine);
	}

	@Test
	public void test() throws Exception {
		depiction.setMolecule("CC(=O)Oc1ccncc1C(O)=O");
		depiction.setDims(800, 800);
		depiction.setHideHydrogens(true);
		depiction.setAromatize(true);
		if (show) {
			showImage(depiction.getImage());
		}
	}
}
