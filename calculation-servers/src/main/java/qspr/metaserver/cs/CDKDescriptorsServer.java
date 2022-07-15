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

package qspr.metaserver.cs;

import com.eadmet.utils.OCHEMUtils;

import qspr.dao.ChemDAO;
import qspr.dao.ChemInfEngine;
import qspr.dao.Various;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsCDKConfiguration;
import qspr.workflow.datatypes.DataTable;

abstract public class CDKDescriptorsServer extends DescriptorsAbstractExecutableServer {

	public String [] jars;

	@Override
	protected
	DataTable calculateDescriptors(DataTable dtMolecules,
			DescriptorsAbstractConfiguration providedConfiguration, int start,
			int batchSize) throws Exception {

		if(jars == null) {
			if(Various.defaultEngine ==  ChemInfEngine.CDK) { 
					jars = OCHEMUtils.loadJars(javaClassPath);
			}
			else {
				setStatus("Using libraries from: " + getExeFile());
				jars = OCHEMUtils.loadJars(getExeFile());
			}
		}

		AbstractCDKCalculator calc = getWorkingHorse();

		DescriptorsCDKConfiguration cdk = (DescriptorsCDKConfiguration) providedConfiguration;

		if(cdk.getCDKAromatization() != null)
			ChemDAO.aromatizeMolecules(dtMolecules, cdk.getCDKAromatization(), start, batchSize);

		calc.calDescriptors(this, dtMolecules, start, batchSize, cdk);
		DataTable tab = getResults();

		for(int i = 0; i < tab.getRowsSize() ; i++) 
			if(!tab.getRow(i).isError())
			{
				int ok = 0;
				for(int j = 0 ; j< tab.getColumnsSize(); j++)
					if(((Double)tab.getValue(i, j)).isNaN())tab.setValue(i, j, 0.);
					else
						ok++; // do we have at least one non zero descriptor?!
				if(ok < tab.getColumnsSize()/2)tab.getRow(i).setError("Descriptor calculation failed for this molecule.");
			}

		return tab;
	}

	AbstractCDKCalculator getWorkingHorse() throws Exception{
		Class<?> clazz;
		if(Various.defaultEngine ==  ChemInfEngine.CDK) { 
			setStatus("Using CDK2");
			clazz = Class.forName("qspr.metaserver.cs.CDK2Calculator"); // There is no way to use old CDK
		}
		else {
			setStatus("Using old CDK");
			clazz = Class.forName("qspr.metaserver.cs.CDKCalculator"); // Using old CDK
			//N.B.! It will always fail in Eclipse since newer CDK libraries will be initialised there
		}

		return (AbstractCDKCalculator)clazz.newInstance();
	}

	@Override
	int getBatchSize() {
		return 1;
	}


}
