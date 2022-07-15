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

import java.util.List;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.Molecule;

@DatabaseMaintenanceJob
public class CalulateInchiTask extends OchemCronjobTask {

	@SuppressWarnings("unchecked")
	public void executeTask() throws Exception
	{
		Globals.startMainTransaction();
		List<Molecule> molList = Globals.session().createCriteria(Molecule.class).add(Restrictions.isNull("mapping1")).list();
		if(molList.size() > 0)
		{
			for (Molecule molecule : molList) 
			{
				molecule.updateInchi();
				Globals.session().saveOrUpdate(molecule);
			}
		}
		Globals.commitMainTransaction();
	}

	public static void main(String[] args)
	{
		new CalulateInchiTask().executeInternal();
	}
}
