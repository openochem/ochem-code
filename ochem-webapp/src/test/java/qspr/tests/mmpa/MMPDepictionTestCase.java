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

package qspr.tests.mmpa;

import org.junit.Test;

import com.eadmet.mmpa.MMPDepiction;

import qspr.dao.ChemInfEngine;
import qspr.tests.depiction.DepictionTestCaseBase;

public class MMPDepictionTestCase extends DepictionTestCaseBase {

	public MMPDepictionTestCase(ChemInfEngine engine) throws Exception {
		super(engine);
		// TODO Auto-generated constructor stub
	}

	@Test
	public void test() throws Exception {
		depiction.setDims(100, 100);
		depiction.setAlpha(0.5);
		depiction.setFormat("png");
		depiction.setColor("FFFFFF");
		depiction.setErrorColor("F5A9A9");
		depiction.setHideHydrogens(false);
		depiction.setAromatize(true);
		depiction.setMolecule("CC(=O)Oc1ccccc1C(O)=O");
		
		if (show) {
			showImage(depiction.getImage());
		}
		
		MMPDepiction depictionMMP = MMPDepiction.get(
				depictionChoice,
				"CC(=O)Oc1cc(CCN(C)C)ccc1C(O)=O",
				"CC(=O)Oc1ccccc1C(O)=O"
		);
		depictionMMP.configure(depiction);
		depictionMMP.search();
		depictionMMP.setDims(300, 300);
		byte[] imgBinary = depictionMMP.getImage();
		
		if (show) {
			showImage(imgBinary);
		}
	}

}
