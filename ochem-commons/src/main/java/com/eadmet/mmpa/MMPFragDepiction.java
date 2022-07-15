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

package com.eadmet.mmpa;

import java.io.IOException;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.dao.ChemDAO;
import qspr.dao.ChemInfEngine;
import qspr.depiction.MoleculeDepiction;

public abstract class MMPFragDepiction extends MoleculeDepiction {

	static MMPFragDepiction cdk,chem;

	public static MMPFragDepiction get(ChemInfEngine impl) {
		try {
			switch (impl) {
			case CDK:
				if(cdk == null)cdk = (MMPFragDepiction) Class.forName("com.eadmet.mmpa.MMPFragDepictionCDK").newInstance();
				return cdk;

			case CHEMAXON:
				if(chem == null)chem = (MMPFragDepiction) Class.forName("com.eadmet.mmpa.MMPFragDepictionChemAxon").newInstance();
				return chem;

			default:
				throw new UserFriendlyException("MMP Depiction implementation unavailable: " + impl);
			}
		} catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
			throw new UserFriendlyException(e);
		}
	}

	//	public static MMPFragDepiction get() {
	//		return get(defaultImpl);
	//	}

	//	public byte[] getImage(MMPFragment frag) throws IOException {
	//		setMolecule(frag.smiles);
	//		return getImage();
	//	}

	public byte[] getImage(String smiles) throws IOException {
		try {
			setMolecule(smiles);
			return getImage();
		}catch(Throwable e) {
			for(ChemInfEngine engine :ChemInfEngine.values())
				if(!ChemDAO.ignoreEngine(engine)) {
					MMPFragDepiction depiction = get(engine);
					depiction.setMolecule(smiles);
					return depiction.getImage();
				}
			throw e;
		}
	}
}
