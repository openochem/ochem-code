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

import org.hibernate.Query;
import org.hibernate.criterion.Restrictions;

import com.eadmet.utils.mailer.Mailer;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.entities.User;

@DatabaseMaintenanceJob
public class CCBYLicenseUpdateTask extends OchemCronjobTask {


	@SuppressWarnings("unchecked")
	public void executeTask() throws Exception
	{
		Globals.startMainTransaction();
		List<User> userList = Globals.session().createCriteria(User.getCurrentClass()).add(Restrictions.eq("licenseVersion",2)).list();
		if(userList.size() > 0)
		{
			for (User user : userList)try
			{
				Query query = Globals.session().createSQLQuery("update ExperimentalProperty set rights=3 where rights=2 and introducer_id="+user.id);
				int records = query.executeUpdate();
				user.licenseVersion = 3;
				Globals.session().saveOrUpdate(user);
				if(records > 0 && !OCHEMConfiguration.inhouseInstallation)
					Mailer.notifyDevelopers("OCHEM CC-BY update", "User " + user.login +" has updated license and " + records + " were updated");
			}
			catch(Exception e){
				if(!OCHEMConfiguration.inhouseInstallation)
					Mailer.notifyDevelopers("OCHEM CC-BY update failure", "User update " + user.login + "\n\n" + e.getMessage());
			}
		}
		Globals.commitMainTransaction();
	}

	public static void main(String[] args)
	{
		new CCBYLicenseUpdateTask().executeInternal();
	}
}
