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

import java.io.IOException;

import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.smarts.SmartsPattern;

import qspr.dao.Various;
import qspr.util.CDKUtils;
import qspr.workflow.utils.QSPRConstants;

public class ExtendedSMARTCDK extends ExtendedSMART {

	public ExtendedSMARTCDK(String exSmart) {
		super(exSmart);
	}
	
	private static IAtomContainer prepare(String mol) throws IOException {
		String sdf = Various.molecule.convertToFormat(mol, QSPRConstants.SDFAROM_GENERAL_WITHH);
		IAtomContainer ret = CDKUtils.readOneMoleculeInAnyFormat(sdf);
		return ret;
	}

	@Override
	protected int getMatchCount(String molecule, String smart) throws Exception {
		IAtomContainer mol = prepare(molecule);
		SmartsPattern pattern = SmartsPattern.create(smart);
		pattern.setPrepare(false);
		int matches = pattern.matchAll(mol).countUnique();
		return matches;
//		int[] matches = pattern.match(mol);
//		return matches.length;
	}

	@Override
	public boolean match(String molecule, String smart) throws Exception {
		if (getMatchCount(molecule, smart) > 0) {
			return true;
		} else {
			return false;
		}
	}

}
