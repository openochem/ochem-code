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

package qspr.metaserver.cs.util;

import java.io.IOException;
import qspr.dao.Various;
import qspr.metaserver.util.MoleculeStandartizer;
import qspr.workflow.utils.SDFProcessor;

import com.eadmet.utils.OCHEMUtils;

public class SavedMolecule
{
	public String md5;
	public int mapping2;    

	public String externId;
	public String originalStructure;
	public String processedStructure;
	public String optimizedStructure;
	public String error;
	public boolean cached = false;
	static MoleculeStandartizer standard;

	public SavedMolecule(String sdf, String externalId, int id) throws IOException
	{
		this.externId = externalId;
		if(sdf == null) return;  // can be null!
		originalStructure = SDFProcessor.standartize(sdf);
		int size = originalStructure.length();
		originalStructure = originalStructure.replaceAll("-0.0000", " 0.0000");
		if(size != originalStructure.length())throw new IOException("Failed elimination of 0 for " + id);
		processedStructure = Various.molecule.addHydrogensAndRemoveRadicalsAndSMARTS(sdf);
		md5 = OCHEMUtils.getMD5(originalStructure);
		mapping2 = id;
	}

	private String removeNH3(String smiles) throws IOException {
		smiles = Various.molecule.convertToCanonicalSMILES(smiles);
		smiles = smiles.replaceAll("\\[NH3", "[N");
		smiles = smiles.replaceAll("\\[NH2", "[N");
		return smiles;
	}

	public String compareFormulas() {
		try {
			if(standard==null)standard = MoleculeStandartizer.getInstance();
			standard.setStandardize();
			return Various.molecule.compareMolecules(removeNH3(processedStructure), removeNH3(standard.doStandartization(optimizedStructure)));
		}catch(Exception e) {
			return "Comparison by formula failed: " + e.getMessage();
		}
	}

	public static void main(String[] args) throws Exception {
		SavedMolecule m = new SavedMolecule("O=C1CCCCC[NH2+]1","",0);
		m.optimizedStructure="O=C1CCCCC[N+]1";
		System.out.println(m.compareFormulas());
	}
}
