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

package com.eadmet.business;

import org.apache.commons.lang.NotImplementedException;

import com.eadmet.exceptions.UserFriendlyException;

//import qspr.OCHEMConfiguration;
import qspr.dao.ChemInfEngine;
//import qspr.entities.AlertSubstitutionVariable;
import qspr.metaserver.util.ExtendedSMART;

public interface SmartsMatcher {
	
	public boolean match(String smiles, String smarts);
	
	public static SmartsMatcher get(ChemInfEngine engine) {
		try {
			switch (engine) {
			case CHEMAXON:
				return (SmartsMatcher) Class.forName("com.eadmet.business.SMARTSMatcherChemaxon").newInstance();
			case CDK:
				return (SmartsMatcher) Class.forName("com.eadmet.business.SMARTSMatcherCDK").newInstance();
			default:
				throw new NotImplementedException("Unknown engine:" + engine.toString());
			}
		} catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
			throw new UserFriendlyException(e);
		}
		
	};
	
//	public static SmartsMatcher get() {
//		return get(OCHEMConfiguration.getCheminfEngine());
//	}
	
	public static ExtendedSMART getExtended(String smarts, ChemInfEngine engine) {
		return ExtendedSMART.create(smarts, engine);
	}
	
	public static ExtendedSMART getExtended(String smarts, ChemInfEngine engine, boolean check) {
		ExtendedSMART exs = SmartsMatcher.getExtended(smarts, engine);
		if (exs.invalid)
			throw new UserFriendlyException("Invalid SMART pattern!");
		return exs;
	}
}
