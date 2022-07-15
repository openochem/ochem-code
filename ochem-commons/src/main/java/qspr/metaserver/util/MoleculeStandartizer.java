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

import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.eadmet.exceptions.CriticalException;

import qspr.dao.ChemInfEngine;
import qspr.dao.Various;
import qspr.metaserver.CalculationServer;
import qspr.metaserver.Event;
import qspr.workflow.datatypes.DataTable;

// Midnighter
// An abstract stub-class for molecule standartization
public abstract class MoleculeStandartizer
{
	boolean deSalt = false, deAromatize = false, neutralizeIt = false, standardizeIt = false, cleanIt = false, deConvertNO2=false;
	public boolean addExplicitHydrogens = false;
	public MolFormat outFormat = MolFormat.SDF;

	public Event<String> statusChange = new Event<String>(this);

	public void setStatus(String status)
	{
		logger.info(status);
		statusChange.fire(status);
	}

	public void setDesalt()
	{
		deSalt = true;
	}

	public void setDearomatize()
	{
		deAromatize = true;
	}

	public void setNeutralize()
	{
		neutralizeIt = true;
	}

	public void setStandardize()
	{
		standardizeIt = true;
	}

	public void setCleanStructure()
	{
		cleanIt = true;
	}

	public void setDeConvertNO2()
	{
		deConvertNO2 = true;
	}

	public static MoleculeStandartizer getInstance(ChemInfEngine engine, String serverHomeDir) throws Exception
	{	
		if(Various.molecule.engine != engine) // at this moment Various.molecule.engine should be available
			throw new CriticalException("Various.molecule.engine is:" + Various.molecule.engine + ", which is different from supplied engine: " + engine.name());

		switch (engine) {
		case CDK:
			logger.info("Using CDK version -> " + Class.forName("org.openscience.cdk.exception.CDKException").getProtectionDomain().getCodeSource().getLocation()
					.toURI());
			return (MoleculeStandartizer) Class.forName("qspr.metaserver.util.CDKStandartizer").newInstance();
		case CHEMAXON:
			return (MoleculeStandartizer) Class.forName("qspr.metaserver.util.ChemaxonStandartizer").newInstance();
		default:
			throw new RuntimeException("Unsupported standartizer: " + engine.name());
		}
	}
	
	public static MoleculeStandartizer getInstance() throws Exception {
		return getInstance(Various.molecule.engine(), null);
	}

	public String doStandartization(String sdf) throws Exception{
		List<String> a = new ArrayList<String>();
		a.add(sdf);
		a = doStandartization(a);
		return a.get(0);
	}

	abstract public List<String> doStandartization(List<String> sdfs) throws Exception;

	abstract public DataTable doStandartizationTable(DataTable datatable, CalculationServer server);

	private static final Logger logger = LogManager.getLogger(MoleculeStandartizer.class);

	public enum MolFormat {SDF, SMILES};  
}

