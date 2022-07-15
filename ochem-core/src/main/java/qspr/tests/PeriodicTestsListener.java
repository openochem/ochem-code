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

package qspr.tests;

import java.io.File;
import java.io.FileWriter;
import java.lang.reflect.Method;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.runner.Description;
import org.junit.runner.Result;
import org.junit.runner.notification.Failure;
import org.junit.runner.notification.RunListener;

import qspr.Globals;
import qspr.OCHEMConfiguration;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.MAILERConstants;
import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;

public class PeriodicTestsListener extends RunListener {

	protected Map<String, PeriodicTestResult> runningTests;
	//public Logger logger;
	public boolean saveToDB;
	public String currentlyRunningTest;
	private int finishedTests = 0;

	public PeriodicTestsListener(boolean saveToDB)
	{
		logger.info("saveToDB is " + saveToDB);
		this.saveToDB = saveToDB;
	}

	public PeriodicTestsListener(String testType, boolean saveToDB)
	{
		logger.info("saveToDB is " + saveToDB);
		this.saveToDB = saveToDB;
	}

	@Override
	public void testRunStarted(Description description) 
	{
		runningTests = new HashMap<String, PeriodicTestResult>();
		logger.info("Test run started");
	};

	@Override
	public void testStarted(Description description) 
	{
		currentlyRunningTest = getTestName(description);
		// FIXME: Define the test type here
		PeriodicTestResult tr = new PeriodicTestResult();
		getTestMetadata(tr, description);
		tr.startTime = new Timestamp(Calendar.getInstance().getTimeInMillis());
		tr.name = getTestName(description);
		tr.className = description.getClassName();
		tr.methodName = description.getMethodName();
		runningTests.put(getTestName(description), tr);
		logger.info(tr.name + " started");
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	private void getTestMetadata(PeriodicTestResult pr, Description description)
	{

		String methodName = description.getMethodName(); //tokenizer.nextToken();
		// method names from parameterized tests look like testMethodName[identifier]
		// to get method annotation the [...] part is removed for identification
		if (methodName.contains("["))
			methodName = methodName.replaceFirst("\\[.+\\]", "");

		String className = description.getClassName(); //tokenizer.nextToken(")").substring(1);
		if (className.startsWith("org.junit"))
			return;

		try
		{
			Class clazz = PeriodicTestRunner.loadClass(className);
			Method m = clazz.getMethod(methodName);

			// Identify the test type
			PeriodicTest methodAnnotation = m.getAnnotation(PeriodicTest.class);
			PeriodicTest classAnnotation = (PeriodicTest) clazz.getAnnotation(PeriodicTest.class);
			if (classAnnotation != null)
				pr.testType = classAnnotation.type();
			if (methodAnnotation != null && methodAnnotation.type() != null && !"".equals(methodAnnotation.type()))
				pr.testType = methodAnnotation.type();
			if (pr.testType == null || pr.testType.equals(""))
				throw new Exception("Unspecified test type for test " + description.getMethodName());

			// Identify the failure notification strategy
			NotifyOnFailure nofMethodAnnotation = m.getAnnotation(NotifyOnFailure.class);
			NotifyOnFailure nofClassAnnotation = (NotifyOnFailure) clazz.getAnnotation(NotifyOnFailure.class);

			if (nofClassAnnotation != null)
				pr.minConsequentFailures = nofClassAnnotation.minConsequentFailures();
			if (nofMethodAnnotation != null)
				pr.minConsequentFailures = nofMethodAnnotation.minConsequentFailures();

		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	public String getTestName(Description description)
	{
		String testName = description.getClassName().replaceAll("qspr.tests.", "").replaceAll("webDriver.", "") + "_" + description.getMethodName();
		return testName;
	}

	@Override
	public void testFailure(Failure failure) 
	{
		PeriodicTestResult tr = runningTests.get(getTestName(failure.getDescription()));
		if (tr == null)
		{
			logger.warn(getTestName(failure.getDescription()) + " was not started properly. The test has been ignored. Failure message: " + failure.getMessage());
			return;
		}

		tr.succeeded = false;

		// If its a user-friendly exception, show a nice message rather than the full stack trace
		if (failure.getException() != null && failure.getException() instanceof UserFriendlyException)
			tr.detailedStatus = failure.getMessage();
		else
			tr.detailedStatus = failure.getTrace(); 
		logger.warn(tr.name + " failed (" + tr.detailedStatus + ")");
		logger.warn(getTestName(failure.getDescription()) + " FAILED: ");
		logger.warn(failure.getMessage());
	}

	@Override
	public void testFinished(Description description) 
	{
		PeriodicTestResult tr = runningTests.get(getTestName(description));
		if (tr == null)
		{
			logger.warn(getTestName(description) + " was not started properly. The test has been ignored.");
			return;
		}
		tr.finishTime = new Timestamp(Calendar.getInstance().getTimeInMillis());
		logger.info(getTestName(description) + " finished");
		logger.info(++finishedTests + " tests are finished.");
	}

	@Override
	public void testRunFinished(Result result) {
		try {
			connectDB();
		} catch (Exception e1) {
			e1.printStackTrace();
		}

		List<String> failures = new ArrayList<String>();
		for (PeriodicTestResult tr : runningTests.values()) {
			if (saveToDB)
				tr.saveToDb(conn);
			if (!tr.succeeded)
				failures.add(tr.name);
		}

		sendEmailNotifications();
		closeConnection();


		if (OCHEMConfiguration.testing) {
			// Send an email with a global test report
			String msg = "Dear OCHEM developer,\n\nA full test run on " + OCHEMConfiguration.rootHost + " has finished.\n";

			if (!failures.isEmpty())
				msg += "\nIn total, <b>"+failures.size()+" tests</b> out of "+runningTests.values().size()+" have failed:\n" + StringUtils.join(failures, "\n");
			else
				msg += "No tests out of " + runningTests.values().size() + " have failed!\n";

			// Temporarily hardcode the emails
			Mailer.postMailSafely(new Email(MAILERConstants.EMAIL_ADMIN, "Full test run finished: " + result.getFailureCount() + " tests failed", msg).useHTML());
		}

		logger.info( (runningTests.isEmpty() ? 0 : runningTests.size() - result.getFailureCount()) + " tests have finished successfully, " + failures.size() + " tests failed");

	}

	@SuppressWarnings("unchecked")
	private void sendEmailNotifications()
	{
		// Send email notifications, if necessary
		Globals.startAllTransactions();
		try
		{
			for (PeriodicTestResult tr : runningTests.values()) 
			{
				if (!tr.succeeded && tr.minConsequentFailures > 0)
				{
					boolean sendNotification = true;
					if (tr.minConsequentFailures > 1)
					{
						// Get the number of the latest consequent failures from DB
						List<Byte> succeededFlags = Globals.session().createSQLQuery("select succeeded from PeriodicTestResult where name=:name order by ptr_id desc limit " + tr.minConsequentFailures).setParameter("name", tr.name).list();
						int failures = 0;
						for (Byte succeededFlag : succeededFlags)
							failures += (1 - succeededFlag);
						sendNotification = (failures >= tr.minConsequentFailures);
					}

					if (sendNotification)
					{
						String suffix = tr.minConsequentFailures > 1 ? " at least " + tr.minConsequentFailures + " times in a row" : "";
						Mailer.notifyDevelopers("OCHEM test \"" + tr.name + "\" failed!", "Regretfully, there is some disturning news.\nThe test " + tr.name + " has failed"+suffix+".\nThis might require your kind attention.\n\nThe system status page can be accessed at " + OCHEMConfiguration.getRootHost() + "/systemstatus/show.do");
					}
				}
			}
			Globals.commitAllTransactions();
		}
		catch (Exception e)
		{
			e.printStackTrace();
			Globals.rollbackAllTransactions();
		}
	}

	private int connectAttemps = 0;
	private static Connection conn; 
	private void connectDB() throws SQLException, InstantiationException, IllegalAccessException, ClassNotFoundException
	{
		try 
		{
			Class.forName(Globals.ochemConf.getProperty("hibernate.connection.driver_class")).newInstance();
			conn = DriverManager.getConnection(Globals.ochemConf.getProperty("hibernate.connection.url"), Globals.ochemConf.getProperty("hibernate.connection.username"), Globals.ochemConf.getProperty("hibernate.connection.password"));
		} 
		catch (SQLException e) 
		{
			if (connectAttemps++ < 3)
				closeConnection();
			e.printStackTrace();
			logger.info("Could not connect to DB to save the test result");
		}
	}



	private int closeAttemps = 0;
	private void closeConnection() {
		try 
		{
			conn.close();
		} 
		catch (SQLException e) 
		{
			if (closeAttemps++ < 3)
				closeConnection();
			e.printStackTrace();
			logger.info("Could not connect to DB to save the test result");
		}
	}

	private static Logger logger = LogManager.getLogger(PeriodicTestsListener.class);
}

class PeriodicTestResult 
{
	public String name;
	public boolean succeeded = true;
	public String className;
	public String methodName;
	public String detailedStatus;
	public Timestamp startTime;
	public Timestamp finishTime;
	public String testType;
	public int minConsequentFailures = 0;

	public void saveToDb(Connection conn) {
		try {
			PreparedStatement ps = conn.prepareStatement("insert into PeriodicTestResult(ptr_id, name, succeeded, detailed_status, start_time, finish_time, test_type, class_name, method_name) values(?, ?, ?, ?, ?, ?, ?, ?, ?)");
			ps.setInt(1, 0);
			ps.setString(2, name);
			ps.setBoolean(3, succeeded);
			ps.setString(4, detailedStatus);
			ps.setTimestamp(5, startTime);
			ps.setTimestamp(6, finishTime);
			ps.setString(7, ArrayUtils.toString(testType));
			ps.setString(8, className);
			ps.setString(9, methodName);
			ps.execute();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void saveToText(String path) {
		try {
			File f = new File(path + "/" + name + "-" + startTime + ".log");
			FileWriter fw = new FileWriter(f);
			fw.write(name + "\n");
			fw.write(succeeded + "\n");
			fw.write(detailedStatus + "\n");
			fw.write(startTime + "\n");
			fw.write(finishTime + "");
			fw.flush();
			fw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}



}