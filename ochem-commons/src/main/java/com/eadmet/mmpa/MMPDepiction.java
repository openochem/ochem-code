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

import java.awt.Color;
import java.io.IOException;

import com.eadmet.exceptions.UserFriendlyException;

import qspr.dao.ChemInfEngine;
import qspr.depiction.MoleculeDepiction;

public abstract class MMPDepiction extends MoleculeDepiction {
	protected Color mmpColor = new Color(255, 162, 0);
	
	protected MMPDepiction() {
		// no action, hidden
	}
	
	public MMPDepiction(String molecule, String template) {
		initialize(molecule, template);
	}
	
	abstract public void initialize(String molecule, String template);
	abstract public void search();
	
	public static MMPDepiction get(ChemInfEngine choice) {
		try {
			switch (choice) {
		        case CDK:
				try {
					return (MMPDepiction) Class.forName("com.eadmet.mmpa.MMPDepictionCDK").newInstance();
				} catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
		                
		        case CHEMAXON:
		        	return (MMPDepiction) Class.forName("com.eadmet.mmpa.MMPDepictionChemAxon").newInstance();
		                    
		        default:
		            throw new UserFriendlyException("MMP Depiction implementation unavailable: " + choice);
		    }
		} catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
			throw new UserFriendlyException(e);
		}
	};
	
	public static MMPDepiction get(ChemInfEngine choice, String molecule, String template) {
		MMPDepiction impl = get(choice);
		impl.initialize(molecule, template);
		return impl;
	}
	
//	public static MMPDepiction get(String molecule, String template) {
//		return get(defaultImpl, molecule, template);
//	}
	
	public Color getHighlightColor() {
		return mmpColor;
	}

	public void setHighlightColor(Color mmpColor) {
		this.mmpColor = mmpColor;
	}
	
	@Override
	public void setMolecule(String mol) throws IOException {
		initialize(mol, mol);
	}
}
