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

import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.mail.MessagingException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.Transactional;
import qspr.business.PendingTaskPeer;
import qspr.entities.Model;
import qspr.entities.User;
import qspr.util.DynaWrap;

import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;

// Delete the models that have not been accessed for more than 6 months (3 months for pending tasks)
// Warn users two times: 2 weeks and 1 week before deletion
// Midnighter

// TODO: Code should be updated to treat PendingTasks (also prediction tasks) properly (upto 1 day work)

@DatabaseMaintenanceJob
public class CleanModelsTask extends OchemCronjobTask
{

	private static final Logger logger = LogManager.getLogger(CleanModelsTask.class);

	boolean debug = false;
	boolean disableScheduler = false;

	public void executeTask() throws Exception
	{
		if (!OCHEMConfiguration.deleteOldModels)
		{
			log("Deleting of models disabled - skipping");
			return;
		}

		log("Starting task deletion of old models");

		Globals.startMainTransaction();

		GregorianCalendar gCal = new GregorianCalendar();
		gCal.add(Calendar.DATE, -7);
		Timestamp weekAgo = new Timestamp(gCal.getTimeInMillis());
		Timestamp lastRun = (Timestamp) Globals.session().createQuery("select max(deletionWarningSent) from Model").uniqueResult();
		if (lastRun == null || lastRun.before(weekAgo) || debug || disableScheduler)
		{
			processModels(false);
			processModels(true);
			processPendingTasks();
		} else
			log("No need to run, last run was just on " + lastRun);
		Globals.commitMainTransaction();
	}

	@SuppressWarnings("unchecked")
	private void processModels(boolean pending) throws MessagingException
	{
		GregorianCalendar gCal = new GregorianCalendar();
		gCal.add(Calendar.MONTH, -6);
		Timestamp sixMonthsAgo = new Timestamp(gCal.getTimeInMillis());

		gCal = new GregorianCalendar();
		gCal.add(Calendar.MONTH, -3);
		Timestamp threeMonthAgo = new Timestamp(gCal.getTimeInMillis());

		gCal = new GregorianCalendar();
		gCal.add(Calendar.DATE, -7);
		Timestamp weekAgo = new Timestamp(gCal.getTimeInMillis());
		Timestamp now = new Timestamp(Calendar.getInstance().getTimeInMillis());

		// Task with id 0 should not exist 
		Globals.session()
		.createQuery("update Model m set m.taskId = null where m.taskId = 0")
		.executeUpdate();

		String condition = pending ? "m.taskId > 0" : "m.taskId is null";
		List<Object[]> results2weeks = Globals.session()
				.createQuery("select m.publicId, s.user, m.name  from Model m join m.session s where (m.lastAccess <= :longTimeAgo) and (m.deletionWarningSent is null) and (m.published = 0 ) and ("+condition+") order by s.user")
				.setTimestamp("longTimeAgo", pending ? threeMonthAgo : sixMonthsAgo).list();

		List<Object[]> results1week = Globals.session()
				.createQuery("select m.publicId, s.user, m.name from Model m join m.session s where (m.deleteAfterDays = 14) and (m.deletionWarningSent <= :weekAgo) and (m.published = 0 ) and ("+condition+") order by s.user")
				.setTimestamp("weekAgo", weekAgo).list();

		List<Object[]> resultsNow = Globals.session()
				.createQuery("select m.publicId, s.user, m.name from Model m join m.session s where (m.deleteAfterDays = 7) and (m.deletionWarningSent <= :weekAgo) and (m.published = 0 ) and ("+condition+") order by s.user")
				.setTimestamp("weekAgo", weekAgo).list();

		Map<User, List<Model>> models2weeks = getModelsGroupedByUsers(results2weeks);
		Map<User, List<Model>> models1week = getModelsGroupedByUsers(results1week);
		Map<User, List<Model>> modelsNow = getModelsGroupedByUsers(resultsNow);

		Set<User> affectedUsers = new HashSet<User>();
		affectedUsers.addAll(models2weeks.keySet());
		affectedUsers.addAll(models1week.keySet());
		affectedUsers.addAll(modelsNow.keySet());

		log("Cleaning after 6 months:" + sixMonthsAgo +" and 3 months: " + threeMonthAgo + " pending = " + pending + " users: " + affectedUsers.size());

		for (User user : affectedUsers)
		{
			DynaWrap  extended  = (user.isExtended())  ?  user.dynaWrapped() :  null;			
			String email = pending ?
					String.format(
							"Dear user, \n\nThank you for creating models on OCHEM!\n\nWe would like to notify you that some of your models have been trained and are awaiting for your decision (save or discard).\n"+
									"OCHEM stores pending models for one month since the last access. The pending models that have not been saved for more than 1 month will be automatically deleted.\n" + 
									"To prevent the deletion, please access your pending models them on the ochem.eu website:<br/><a class='fb-button' href='"+OCHEMConfiguration.rootHost+"/user/profile.do?login=%s'>Open your OCHEM profile and review your pending models</a>.\n\n"+
									"The saved models will be still kept for next 6 months. \n\n",
									user.login)
									:
										String.format(
												"Dear user, \n\nThank you for creating models on OCHEM!\n\nWe would like to notify you that some of your models have not been accessed for a long time and, therefore, will be automatically deleted.\n"+
														"To prevent the deletion, please access them on the ochem.eu website or choose to keep all your models:<br/><a class='fb-button' href='"+OCHEMConfiguration.rootHost+"/user/profile.do?login=%s'>Open your OCHEM profile and review your models</a>.\n\n"+
														"The accessed models will be still kept for next 6 months. \n\n"+
														"Please note that if you publish your model, it will never expire and will be never automatically deleted.\n\n",
														user.login);

							if (modelsNow.get(user) != null)
								email += "For "+modelsNow.get(user).size()+" model(s), we have sent you a notification two weeks ago, but you did not access these models since that time.\nTherefore, the following models will be deleted NOW: \n" + getModelLinks(modelsNow.get(user)) + "\n\n";

							if (models1week.get(user) != null)
								email += "The following models will be deleted in ONE week:\n" + getModelLinks(models1week.get(user)) + "\n\n";

							if (models2weeks.get(user) != null)
								email += "The following models will be deleted in TWO weeks:\n" + getModelLinks(models2weeks.get(user));

							email += "\n\n_______________________________\nThis e-mail was generated automatically by the OCHEM database. In case of any questions, please reply to this email";

							try
							{
								log(email);

								if (!debug)
								{
									
									if (extended !=  null) {
										try
										{
											Mailer.postMail(new Email(extended.getString("email"), "Your unused models on the OCHEM will be deleted", email).useHTML());
										}
										catch (Exception e)
										{
											// If the users email does not accept our messages, its his problem. Its no reason to fail.
											log("WARNING: Could not send email to " + extended.getString("email"));
										}
									}
									
									// 2 weeks notifications
									if (models2weeks.get(user) != null)
										Globals.session()
										.createQuery("update Model m set m.deletionWarningSent = :now, m.deleteAfterDays = 14 where m.publicId in (:ids)")
										.setParameterList("ids", getModelIdentifiers(models2weeks.get(user)))
										.setTimestamp("now", now)
										.executeUpdate();

									// 1 week notifications
									if (models1week.get(user) != null)
										Globals.session()
										.createQuery("update Model m set m.deletionWarningSent = :now, m.deleteAfterDays = 7 where m.publicId in (:ids)")
										.setParameterList("ids", getModelIdentifiers(models1week.get(user)))
										.setTimestamp("now", now)
										.executeUpdate();

									// Delete the models
									if (modelsNow.get(user) != null)
									{
										List<Long> userModelIDs = getModelIdentifiers(modelsNow.get(user));
										List<Long> mmIds = Globals.session().createQuery("select mmm.id from ModelMapping mmm join mmm.model m where m.publicId in (:ids)")
												.setParameterList("ids", userModelIDs)
												.list();

										List<Long> pTaskIds = Globals.session().createQuery("select pt.id from PendingTask pt join pt.model m where m.publicId in (:ids)")
												.setParameterList("ids", userModelIDs)
												.list();

										List<Long> modelIDs = Globals.session().createQuery("select id from Model where publicId in (:ids)")
												.setParameterList("ids", userModelIDs)
												.list();

										if (!pTaskIds.isEmpty())
											Globals.session()
											.createQuery("delete from PendingTask pt where pt.id in (:ptIds)")
											.setParameterList("ptIds", pTaskIds)
											.executeUpdate();
										if (!mmIds.isEmpty())
											Globals.session()
											.createQuery("delete from ModelMapping mm where mm.id in (:mmIds)")
											.setParameterList("mmIds", mmIds)
											.executeUpdate();
										if (!modelIDs.isEmpty())
										{
											Globals.session().createSQLQuery("delete from ModelTag where model_id in (:modelIDs)")
											.setParameterList("modelIDs", modelIDs)
											.executeUpdate();
											Globals.session().createSQLQuery("delete from ValidationSet where model_id in (:modelIDs)")
											.setParameterList("modelIDs", modelIDs)
											.executeUpdate();
											Globals.session().createSQLQuery("delete from CachedPrediction where model_id in (:modelIDs)")
											.setParameterList("modelIDs", modelIDs)
											.executeUpdate();
										}

										Globals.session()
										.createQuery("delete from Model m where m.publicId in (:ids)")
										.setParameterList("ids", userModelIDs)
										.executeUpdate();
									}
								}

								Globals.session().flush();
							}
							catch (Exception e)
							{
								Mailer.notifyDevelopers(e, "CleanModelTasks");
								e.printStackTrace();
							}
		}
	}

	@SuppressWarnings("unchecked")
	private void processPendingTasks() throws Exception {

		PendingTaskPeer.updateTaskStatuses(null, null, 1000);

		GregorianCalendar gCal = new GregorianCalendar();
		gCal.add(Calendar.MONTH, -3);
		Timestamp threeMonthsAgo = new Timestamp(gCal.getTimeInMillis());

		List<Long> taskIDs = Globals.session().createQuery("select id from PendingTask where type != 'MODEL_TRAINING' and timePosted <= :threeMonthsAgo and published=0")
				.setParameter("threeMonthsAgo", threeMonthsAgo).list();

		logger.info("" + taskIDs.size() + " pending tasks to delete");
		while (!taskIDs.isEmpty()) {
			List<Long> batchIDs = taskIDs.subList(0, Math.min(100, taskIDs.size()));
			logger.info("Deleting " + batchIDs.size() + " pending tasks");
			Globals.session().createQuery("delete from PendingTask where id in (:ids)")
			.setParameterList("ids", batchIDs).executeUpdate();
			batchIDs.clear();
		}

	}

	private String getModelLinks(List<Model> models)
	{
		StringBuilder urls = new StringBuilder();
		for (Model model : models)
			urls.append(OCHEMConfiguration.getRootHost()+OCHEMConfiguration.rootDir + "/model/" + model.publicId + " ("+model.name+")\n");
		return urls.toString();
	}

	private List<Long> getModelIdentifiers(List<Model> models)
	{
		List<Long> ids = new ArrayList<Long>();
		for (Model model : models)
			ids.add(model.publicId);
		return ids;
	}

	// Get a list of models
	private Map<User, List<Model>> getModelsGroupedByUsers(List<Object[]> results)
	{
		Map<User, List<Model>> result = new HashMap<User, List<Model>>();

		User user = null;
		List<Model> models = new ArrayList<Model>();
		for (Object[] objects : results)
		{
			User curUser = (User) objects[1];
			if (!curUser.equals(user))
			{
				if (user != null)
				{
					result.put(user, models);
					models = new ArrayList<Model>();
				}
				user = curUser;
			}

			Model model = new Model();
			model.publicId = (Long) objects[0];
			model.name = (String) objects[2];
			models.add(model);
		}

		if (user != null)
			result.put(user, models);

		return result;
	}

	public static void main(String[] args) throws Exception
	{
		OCHEMConfiguration.mirror = false;

		new Transactional()
		{

			@Override
			protected void wrapped() throws Exception
			{
				Mailer.enable = true;
				new CleanModelsTask().processModels(false);
			}
		}.execute();

		//new CleanModelsTask().postMail("test@gmail.com", "A teeest", "Really.-");
	}
}
