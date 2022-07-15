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
import java.math.BigDecimal;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Map;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Marshaller;
import javax.xml.datatype.DatatypeConfigurationException;
import javax.xml.datatype.DatatypeFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.runner.Description;
import org.junit.runner.Result;
import org.junit.runner.notification.Failure;
import org.junit.runner.notification.RunListener;

import qspr.tests.schema.Testsuite;
import qspr.tests.schema.Testsuite.Testcase;
import qspr.tests.schema.Testsuites;

/**
 * Generates Ant-compatible XML reports for junit tests
 *
 */
public class ReportGeneratorTestListener extends RunListener 
{
	private static Logger logger = LogManager.getLogger(ReportGeneratorTestListener.class);
	private Testsuites suites;
	private File file;
	private Map<Object, Long> startingTimes = new HashMap<Object, Long>();

	public ReportGeneratorTestListener(File file) throws Exception {
		this.file = file;
		logger.info("Initialized with file " + file);
	}

	private void addTestCasesToSuite(Description description, Testsuite suite) {
		if (description.getChildren().size() == 0) {
			suite.setTests(suite.getTests() + 1);
			Testcase testCase = new Testcase();
			testCase.setClassname(description.getClassName());
			testCase.setName(description.getMethodName());
			suite.getTestcase().add(testCase);
			startingTimes.put(testCase, Calendar.getInstance().getTimeInMillis());
		} else {
			for (Description d : description.getChildren())
				addTestCasesToSuite(d, suite);
		}
	}

	private Testcase findTestCase(Description description) {
		for (Testsuite suite : suites.getTestsuite())
			for (Testcase c : suite.getTestcase())
				if (c.getClassname().equals(description.getClassName()) && c.getName().equals(description.getMethodName()))
					return c;
		return null;
	}

	private Testsuite findTestSuite(Description description) {
		for (Testsuite suite : suites.getTestsuite())
			for (Testcase c : suite.getTestcase())
				if (c.getClassname().equals(description.getClassName()) && c.getName().equals(description.getMethodName()))
					return suite;
		return null;
	}

	@Override
	public void testRunStarted(Description description) {
		suites = new Testsuites();
		for (Description d : description.getChildren()) {
			Testsuites.Testsuite suite = new Testsuites.Testsuite();
			try {
				suite.setTimestamp(DatatypeFactory.newInstance().newXMLGregorianCalendar((GregorianCalendar) GregorianCalendar.getInstance()));
			} catch (DatatypeConfigurationException e) {
				e.printStackTrace();
			}
			try {
				suite.setHostname(InetAddress.getLocalHost().getHostName());
			} catch (UnknownHostException e1) {
				suite.setHostname("TestingServer");
			}

			suite.setName(d.getClassName());
			startingTimes.put(suite, Calendar.getInstance().getTimeInMillis());

			for (Description id : d.getChildren())
				addTestCasesToSuite(id, suite);
			suites.getTestsuite().add(suite);
		}
	}

	@Override
	public void testRunFinished(Result result) {
		for (Testsuite suite : suites.getTestsuite()) {
			suite.setSystemErr("");
			suite.setSystemOut("");
			if (startingTimes.get(suite) != null)
				suite.setTime(getRunningTime(suite));
		}
		try {
			JAXBContext context = JAXBContext.newInstance(Testsuite.class, Testsuites.class);
			Marshaller marshaller = context.createMarshaller();
			marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
			marshaller.marshal(suites, file);
		} catch (Exception e) {
			e.printStackTrace();
		}
		logger.info("Report generated successfully");
	}

	@Override
	public void testStarted(Description description) {
		Testsuite s = findTestSuite(description);
		if (startingTimes.get(s) == null)
			startingTimes.put(s, Calendar.getInstance().getTimeInMillis());

		Testcase c = findTestCase(description);
		startingTimes.put(c, Calendar.getInstance().getTimeInMillis());
	}

	@Override
	public void testFinished(Description description) {
		Testcase c = findTestCase(description);
		c.setTime(getRunningTime(c));
	}
	
	private void doFail(Failure failure, Description description) {
		Testsuite s = findTestSuite(description);
		s.setFailures(s.getFailures() + 1);

		Testcase c = findTestCase(description);
		Testcase.Failure f = new Testcase.Failure();
		f.setMessage(failure.getMessage());
		f.setType(failure.getException().getClass().getName());
		f.setValue(failure.getTrace());
		c.setFailure(f);
		c.setTime(getRunningTime(c));
	}
	
	private BigDecimal getRunningTime(Object o)
	{
		return new BigDecimal((Calendar.getInstance().getTimeInMillis() - startingTimes.get(o)) / 1000D);
	}

	@Override
	public void testFailure(Failure failure) {
		if (failure.getDescription().isSuite()) {
		{
			for (Description d : failure.getDescription().getChildren())
				doFail(failure, d);
		}
		} else {
			doFail(failure, failure.getDescription());
		}
	}

	@Override
	public void testIgnored(Description description) {

	}
}