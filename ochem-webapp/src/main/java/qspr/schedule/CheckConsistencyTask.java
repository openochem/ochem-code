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

import javax.mail.MessagingException;

import org.quartz.JobExecutionException;

import qspr.Globals;

import com.eadmet.utils.mailer.Mailer;

// Delete the models that have not been accessed for more than 6 months (1 month for pending tasks)
// Warn users two times: 2 weeks and 1 week before deletion
// Midnighter

// TODO: Code should be updated to treat PendingTasks (also prediction tasks) properly (upto 1 day work)

@DatabaseMaintenanceJob
public class CheckConsistencyTask extends OchemCronjobTask
{

	boolean debug = false;

	public void executeTask() throws Exception
	{
		Globals.startMainTransaction();

		StringBuffer message = new StringBuffer();
		findInvalidModels(message);
		if (message.length() > 0)
			Mailer.notifyDevelopers("OCHEM database inconsistencies", "There is something you should know.\nNamely, there are some inconsistencies in the OCHEM database:\n\n" + message.toString());
		else
			log("No inconsistencies found");

		Globals.commitMainTransaction();
	}

	@SuppressWarnings("unchecked")
	private String findInvalidModels(StringBuffer message)
	{
		List<Object[]> disappearedModels = Globals.session().createSQLQuery("select login, m.name, model_id from Model m left join PendingTask using (model_id) left join Session on (m.session_id=Session.session_id) natural left join User where m.task_id is not null and PendingTask.task_id is null").list();
		if (!disappearedModels.isEmpty())
		{
			message.append("Unfinished models with lost pending tasks:\n\n");
			for (Object[] row : disappearedModels)
				message.append("model " + row[1] + " by " + row[0] + ", model ID " + row[2] + "\n");
			message.append("\n\n");
		}

		List<Object[]> modelsWithoutProperties = Globals.session().createSQLQuery("select name, model_id from Model natural left join ModelMapping where property_id is null and model_template_id != 31").list();
		if (!modelsWithoutProperties.isEmpty())
		{
			message.append("Invalid models without properties:\n\n");
			for (Object[] row : modelsWithoutProperties)
				message.append("model " + row[0] + ", model ID " + row[1] + "\n");
			message.append("\n\n");
		}

		return message.toString();
	}

	public static void main(String[] args) throws JobExecutionException, MessagingException
	{
		new CheckConsistencyTask().executeInternal();
	}
}
