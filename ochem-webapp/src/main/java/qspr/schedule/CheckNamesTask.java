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

package qspr.schedule;

import java.util.ArrayList;
import java.util.List;

import org.hibernate.Query;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.entities.MoleculeName;
import qspr.entities.Property;
import qspr.entities.ValidatedFact;
import qspr.util.NCBI_Utility;

@DatabaseMaintenanceJob
public class CheckNamesTask  extends OchemCronjobTask { 

	static final int NAMESIZE = 1000;

	@Override
	public boolean shouldRun()
	{
		return super.shouldRun() && OCHEMConfiguration.allowExternalServices;
	}

	// create table, e.g. tmp_aa with mapping2_id to be updated
	// delete from  ValidatedFact all facts with this mapping2_id to re-init validation

	@SuppressWarnings("unchecked")
	public void executeTask() throws Exception
	{
		List<String> namesToCheck = new ArrayList<String>();

		Globals.startMainTransaction();

		Query q = 
				Globals.session().createSQLQuery(
						"select distinct mn.* from MoleculeName mn " +
								"join ExperimentalPropertyName epn using (molecule_name_id) " +
								"join ExperimentalProperty ep using (exp_property_id) " +
								"left join ValidatedFact vf on (mn.molecule_name_id = vf.molecule_name_id and vf.validation = "+ValidatedFact.VALIDATED+") " +
								"where vf.validatedfact_id is null " +
								"and ep.property_id !=" +Property.getDummyId())  
				.addEntity(MoleculeName.class)
				.setMaxResults(NAMESIZE);

		List<MoleculeName> nameList = q.list();

		for (MoleculeName molname : nameList) 
			namesToCheck.add(molname.name);

		log("names to check are:\t" + nameList);

		Globals.commitMainTransaction();

		for (String name : namesToCheck)
		{
			Globals.startAllTransactions();
			try
			{
				log("Search for " + name);
				NCBI_Utility.getMoleculeByName(name);
				Globals.commitAllTransactions();
			} catch (Exception e)
			{
				log("Exception in PubChem search");
				e.printStackTrace();
				Globals.rollbackAllTransactions();
			}
		}
	}


	public static void main(String[] args) throws Exception
	{
		new CheckNamesTask().executeTask();
	}


}
