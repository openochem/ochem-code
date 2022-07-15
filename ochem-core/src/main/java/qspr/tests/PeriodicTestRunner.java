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
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.runner.JUnitCore;
import org.junit.runner.Request;
import org.reflections.Reflections;

import qspr.Environment;
import qspr.OCHEMConfiguration;

import com.eadmet.utils.config.ReflectionsSingleton;

/** The universal runner of modular jUnit tests, which can:
 *  - automatically find the annotated testing classes
 *  - determine the necessity to run based on scheduling options
 *  - log the test results into database (see PeriodicTestListener.java)
 *  - notify developers about critical failures (using a notification strategy, see the @NotifyOnFailure annotation and PeriodicTestListener.java)
 * 
 *  @author Midnighter
 */

public class PeriodicTestRunner implements Runnable
{
	private static transient final Logger logger = LogManager.getLogger(PeriodicTestRunner.class);

	public boolean debug;
	public PeriodicTestsListener listener;
	public PeriodicTestFilter filter;
	public static ClassLoader loader;
	private static volatile PeriodicTestRunner instance;
	public boolean saveReport;
	public String reportFileName;


	public PeriodicTestRunner()
	{
		saveReport = OCHEMConfiguration.testing;
		logger.info("Test runner created");
		listener = new PeriodicTestsListener(null, !debug);
		filter = new PeriodicTestFilter(listener);
	}

	public static PeriodicTestRunner getRunningInstance()
	{
		return instance;
	}

	public void run()
	{
		if (instance != null)
		{
			logger.warn("Tests are already running. Concurrent runs are now allowed. Exiting.");
			return;
		}
		if (reportFileName == null)
			reportFileName = Environment.getWebInfDir()+"/TEST-junit-report.xml";
		try
		{
			if (annotatedTestClasses.length == 0)
				return;
			instance = this;
			logger.info("Starting tests");
			JUnitCore jUnit = new JUnitCore();
			jUnit.addListener(listener);
			Request request = Request.classes(annotatedTestClasses).filterWith(filter);
			int testsSize = request.getRunner().testCount(); // This is necessary to invoke the filter and to know whether there are any runnable tests

			if (saveReport)
			{
				try {
					jUnit.addListener(new ReportGeneratorTestListener(new File(reportFileName)));
				} catch (Exception e) {
					e.printStackTrace();
				}
			}

			if (filter.hasTests)
			{
				logger.info("Starting global test run");
				logger.info(testsSize + " tests are registered to run");
				jUnit.run(request);
				logger.info("Finished global test run");
			}
			else
				logger.info("No tests to run");
		}
		catch (Exception e)
		{
			logger.info(e);
		}
		finally
		{
			instance = null;
		}
	}

	@SuppressWarnings("rawtypes")
	public static Class loadClass(String name) throws ClassNotFoundException
	{
		if (loader == null)
			return Class.forName(name);
		else
			return loader.loadClass(name);
	}

	@SuppressWarnings("rawtypes")
	public static Class[] annotatedTestClasses;
	static
	{
		Reflections reflections = ReflectionsSingleton.get();
		Set<Class<?>> annotated = reflections.getTypesAnnotatedWith(PeriodicTest.class);
		annotatedTestClasses = new Class[annotated.size()];
		int i = 0;
		for (Class<?> clazz : annotated) {
			annotatedTestClasses[i++] = clazz;
		}

		System.out.print("Found " + i + " registered testing classes: ");
		for (i = 0; i < annotated.size(); i++)
			System.out.print(annotatedTestClasses[i].getSimpleName() + ", ");
		logger.info("\n");
	}

	public static void main(String[] args) {
		new PeriodicTestRunner().run();
	}
}
