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

import java.lang.reflect.Method;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.runner.Description;
import org.junit.runner.manipulation.Filter;

import qspr.Globals;
import qspr.OCHEMConfiguration;

import com.eadmet.utils.mailer.Mailer;

/**
 *
 * JUnit filter implementation, which:
 *  - connects to the database to identify the last test runs
 *  - selects the tests that must be executed now based on their scheduling options (see the @PeriodicTest annotation)
 *
 * @author Midnighter
 *
 */

public class PeriodicTestFilter extends Filter
{
	private static transient final Logger logger = LogManager.getLogger(PeriodicTestFilter.class);

	PeriodicTestsListener listener;
	public boolean startupRun;
	public boolean hasTests = false;

	/**
	 * Used for skipping the scheduling info, mostly for debugging purposes
	 */
	public boolean disableScheduler = false;

	/**
	 * Run only the specific test type (e.g., general or selenium tests)
	 */
	public String runSpecificTestType; //= "general"; // TODO: Remove "general" when ready with testing

	/**
	 * In case we want to force particular test running
	 */
	public String runSpecificTest;

	public PeriodicTestFilter(PeriodicTestsListener listener)
	{
		this.listener = listener;
	}

	@Override
	public String describe()
	{
		return "Run tests enabled by configuration";
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	@Override
	public boolean shouldRun(Description description)
	{
		String methodName = description.getMethodName();
		if ("null".equals(methodName))
			methodName = null;

		String className = description.getClassName();

		if (description.getChildren().size() > 0)
			return true;

		if (!description.getDisplayName().contains("("))
			return false;

		String testName = methodName;

		Class clazz;
		PeriodicTest classAnnotation;
		PeriodicTest methodAnnotation = null;
		ScheduleType schedule = null;
		String testType = null;
		if (methodName == null)
			return true;

		try
		{
			clazz = PeriodicTestRunner.loadClass(className);

			// method names from parameterized tests look like testMethodName[identifier]
			// to get method annotation the [...] part is removed for identification
			if (methodName.contains("["))
				methodName = methodName.replaceAll("\\[.+\\]", "");

			if (methodName != null)
			{
				Method m = clazz.getMethod(methodName);
				methodAnnotation = m.getAnnotation(PeriodicTest.class);
			}

			classAnnotation = (PeriodicTest) clazz.getAnnotation(PeriodicTest.class);
			if (classAnnotation != null)
			{
				schedule = classAnnotation.schedule();
				testType = classAnnotation.type();
			}
			if (methodAnnotation != null && methodAnnotation.schedule() != null)
			{

				schedule = methodAnnotation.schedule();
				testType = methodAnnotation.type();
			}
		}
		catch (Exception e)
		{
			e.printStackTrace();
			return false;
		}

		if (schedule == null)
			logger.info("WARNING! No schedule provided for the test " + description.getDisplayName() + ". The test will always run.");

		if ( ! timeToRun(listener.getTestName(description), testType, schedule, getLastRun(listener.getTestName(description))))
			return false;

		try
		{
			if (testName.equals("emailTest"))
				return hasTests = Mailer.enable;
			if (testName.equals("ncbiSearchTest"))
				return hasTests = OCHEMConfiguration.allowExternalServices;
			return hasTests = true;
		} catch (Exception e)
		{
			e.printStackTrace();
		}
		return false;
	}

	private Map<String, Date> lastRuns;
	@SuppressWarnings("unchecked")
	private Date getLastRun(String testName)
	{
		if (lastRuns == null || lastRuns.isEmpty())
		{
			try
			{
				lastRuns = new HashMap<String, Date>();
				logger.info("Getting the times of last test runs...");
				Globals.startAllTransactions();
				List<Object[]> rows = Globals.session().createSQLQuery("select name, max(finish_time) last_run from PeriodicTestResult group by name order by last_run desc").list();
				for (Object[] objects : rows) {
					java.sql.Timestamp date = (java.sql.Timestamp) objects[1];
					if (date != null)
						lastRuns.put("" + objects[0], new Date(date.getTime()));
				}
			}
			finally
			{
				Globals.commitAllTransactions();
			}
		}

		return lastRuns.get(testName);
	}

	private boolean timeToRun(String testName, String testType, ScheduleType schedule, Date lastRun)
	{
		// First consider the specific test types or individual tests
		if (runSpecificTestType != null)
		{
			if (runSpecificTestType.toLowerCase().startsWith("not "))
			{
				if (testType.contains(runSpecificTestType.split(" ")[1]))
					return false;
			}
			else if (!testType.contains(runSpecificTestType))
				return false;
		}

		if (runSpecificTest != null && !runSpecificTest.equalsIgnoreCase(testName) && !testName.toLowerCase().contains(runSpecificTest.toLowerCase()))
			return false;

		// Now consider the scheduler
		if (disableScheduler)
			return true;
		if (schedule == null)
			return true;
		if (lastRun == null)
			return true;
		long timeSinceLastRun = Calendar.getInstance().getTimeInMillis() - lastRun.getTime();
		int requiredInterval = schedule.getMilliseconds();

		logger.debug(testName + " last run was " + timeSinceLastRun + "ms. ago ("+lastRun+"), required interval is " + requiredInterval + "ms.");
		if (schedule.equals(ScheduleType.AT_FIVE_AM))
		{
			Calendar cal = Calendar.getInstance();
			cal.set(Calendar.HOUR_OF_DAY, 0);
			cal.set(Calendar.MINUTE, 0);
			cal.set(Calendar.SECOND, 0);
			cal.set(Calendar.MILLISECOND, 0);

			long midnight = cal.getTimeInMillis();
			long fromMidnightTillNow = Calendar.getInstance().getTimeInMillis() - midnight;
			boolean shouldRunToday = timeSinceLastRun > fromMidnightTillNow;
			boolean isScheduleTime = fromMidnightTillNow > requiredInterval;

			// logger.info(testName + " should run today: " + shouldRunToday + " - is after schedule time: " + isScheduleTime);
			return shouldRunToday && isScheduleTime;
		} 
		else if (requiredInterval > 0)
			return timeSinceLastRun > requiredInterval * 0.95;
			else
				return startupRun;
	}
}
