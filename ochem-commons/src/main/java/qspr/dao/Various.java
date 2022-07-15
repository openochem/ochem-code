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

package qspr.dao;

import java.util.HashMap;
import java.util.Map;

import com.eadmet.exceptions.CriticalException;

import qspr.OCHEMConfiguration;

public class Various {

	public static ChemInfEngine defaultEngine = OCHEMConfiguration.getCheminfEngine();
	public static ChemDAO molecule = getDefaultCheminfImpl();
	public static Map<ChemInfEngine,ChemDAO> cache  = null;

	public static ChemDAO getDefaultCheminfImpl(){
		return getCheminfImpl(defaultEngine);
	}

	// required to upload dynamically libraries from a path
	public static ChemDAO getCheminfImpl(ChemInfEngine engine){
		if(cache == null)cache  = new HashMap<ChemInfEngine,ChemDAO>();
		if(cache.containsKey(engine))return cache.get(engine);

		try {
			ChemDAO imp = null;
			switch (engine) {
			case CDK: // Trying from both possible places or assuming that java is part of the code
				//OCHEMUtils.loadJars(QSPRConstants.CDKLIB);
				//OCHEMUtils.loadJars(QSPRConstants.SOURCE + QSPRConstants.CDK2);
				imp = (ChemDAO) Class.forName("qspr.dao.ChemDAOImplCDK").newInstance();
				break;

			case NONE:
				imp = (ChemDAO) Class.forName("qspr.dao.ChemDAOImplNONE").newInstance();
				break;

			default:
				throw new CriticalException("Engine is unavailable: " + engine);			
			}
			cache.put(engine, imp);
			return imp;
		}catch(Exception e) {
			System.err.println("DAO implementation failed to initalize for engine: " + engine + ". Returning null...");
			e.printStackTrace();
			return null;
		}
	}

}
