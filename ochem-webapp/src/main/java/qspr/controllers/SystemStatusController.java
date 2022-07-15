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

package qspr.controllers;

import java.sql.Timestamp;
import java.util.Calendar;
import java.util.Enumeration;
import java.util.GregorianCalendar;
import java.util.List;

import javax.servlet.ServletOutputStream;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.Criteria;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.entities.PeriodicTestResult;
import qspr.frontend.WebModel;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.transport.CSTransport;
import qspr.schedule.MonitoringTask;
import qspr.tests.PeriodicTestRunner;

import com.eadmet.exceptions.UserFriendlyException;

@Controller
public class SystemStatusController extends ControllerWrapper 
{
	public SystemStatusController()
	{
		sessionRequired = true;
	}

	/**
	 * Show system load
	 * @param req
	 * @param res
	 * @return
	 */
	public ModelAndView monitors(HttpServletRequest req, HttpServletResponse res) {
		if (assertParam("gc"))
			System.gc();
		return new WebModel().setTemplate("system-monitors").getModelAndView();
	}

	public synchronized ModelAndView getMonitorData(HttpServletRequest req, HttpServletResponse res) {
		if (assertParam("gc"))
			System.gc();
		int length = 60*10 / 2;
		if (assertParam("length"))
			length = getIntParam("length") / 2;

		return new ModelAndView("json", "object", new Object[]{getMonitorData("memory", length), getMonitorData("load", length)});
	}

	@SuppressWarnings("rawtypes")
	private List getMonitorData(String metrics, int length) 
	{
		List fullMonitorData = null;
		if ("memory".equals(metrics))
			fullMonitorData = MonitoringTask.memoryMonitor;
		else
			fullMonitorData = MonitoringTask.loadMonitor;
		return fullMonitorData.subList(fullMonitorData.size() - Math.min(length, fullMonitorData.size()), fullMonitorData.size());
	}

	/*
	 * Show results of tests
	 */
	@SuppressWarnings("unchecked")
	public ModelAndView show(HttpServletRequest req, HttpServletResponse res) throws Exception {
		if (!Globals.isSuperUser()) return redirect("home/show.do");
		meta(req, res); // to register IP

		int day = -7;

		String days = req.getParameter("days");
		if(days!=null)day = -Integer.parseInt(days);

		GregorianCalendar gCal = new GregorianCalendar();
		long time2 = gCal.getTimeInMillis();
		gCal.add(Calendar.DATE, day);
		gCal.set(Calendar.HOUR_OF_DAY, 0);
		gCal.set(Calendar.MINUTE, 0);
		gCal.set(Calendar.SECOND, 0);
		long time1 = gCal.getTimeInMillis();

		if (assertParam("rewindDays"))
		{
			long rewind = getLongParam("rewindDays") * 24 * 60 * 60 * 1000;
			time1 -= rewind;
			time2 -= rewind;
		}

		List<PeriodicTestResult> results;

		Criteria criteria = Globals.session().createCriteria(PeriodicTestResult.class)
				.add(Restrictions.between("startTime", new Timestamp(time1), new Timestamp(time2)));

		if (req.getParameter("type") != null && req.getParameter("type").equalsIgnoreCase("selenium"))
			criteria.add(Restrictions.like("testType", "selenium%"));
		else
			criteria.add(Restrictions.or(Restrictions.isNull("testType"), Restrictions.eq("testType", "general")));

		criteria.addOrder(Order.asc("succeeded"));
		criteria.addOrder(Order.asc("startTime"));
		results = criteria.list();

		WebModel wm = new WebModel(WebModel.fromList(results)).setTemplate("system-status");
		if (PeriodicTestRunner.getRunningInstance() != null)
			wm.addParam("running", PeriodicTestRunner.getRunningInstance().listener.currentlyRunningTest);
		return wm.getModelAndView();
	}

	/**
	 * Metaserver
	 * @param req
	 * @param res
	 * @return
	 * @throws Exception
	 */

	public ModelAndView meta(HttpServletRequest req, HttpServletResponse res) throws Exception {
		if (!Globals.isSuperUser()) return redirect("home/show.do");
		OCHEMConfiguration.addTrustedAddress(ThreadScope.resolveRemoteAddr(req));
		CSTransport transport = new CSTransport();
		transport.executeCommand(new Command(Command.CL_REGISTER_ADMIN_IP, ThreadScope.resolveRemoteAddr(req)).sid("Dev Script"));
		return redirect(transport.serverURL);
	}

	@SuppressWarnings("unchecked")
	public ModelAndView headers(HttpServletRequest req, HttpServletResponse res) throws Exception 
	{
		ServletOutputStream os = res.getOutputStream();
		Enumeration<String> headers = req.getHeaderNames();
		while (headers.hasMoreElements())
		{
			String name = headers.nextElement();
			os.println("<b>"+name+"</b>: "+req.getHeader(name)+"<br/>");
		}
		os.close();
		return null;
	}

	public ModelAndView forceRun(HttpServletRequest req, HttpServletResponse res) throws Exception {
		if (PeriodicTestRunner.getRunningInstance() != null)
			throw new UserFriendlyException("The tests are already running. Please, wait.");

		PeriodicTestRunner runner = new PeriodicTestRunner();
		runner.filter.disableScheduler = true;
		if (!"all".equals(getParam("test"))) {
			runner.filter.runSpecificTest = getParam("test");
		}

		if (assertParam("type"))
			runner.filter.runSpecificTestType = getParam("type");

		new Thread(runner).start();

		return redirect("systemstatus/show.do");
	}

}
