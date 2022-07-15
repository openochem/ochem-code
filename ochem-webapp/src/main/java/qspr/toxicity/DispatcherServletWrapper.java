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

package qspr.toxicity;

/*
 	DispatcherServletWrapper
 	Wrapper around Spring's DispatcherServlet
 	Necessary mostly for explicit transaction demarcation 
 */

import java.io.Writer;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.List;

import javax.servlet.ServletConfig;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.servlet.http.HttpSessionEvent;
import javax.servlet.http.HttpSessionListener;

import org.apache.commons.lang.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Criteria;
import org.hibernate.criterion.Restrictions;
import org.springframework.web.servlet.DispatcherServlet;

import qspr.Environment;
import qspr.Globals;
import qspr.LoadLogger;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.entities.Model;
import qspr.entities.Molecule;
import qspr.entities.Session;
import qspr.metaserver.transport.TransportFactory;
import qspr.tests.PeriodicTestRunner;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.MemoryUtils;

@SuppressWarnings("unchecked")
public class DispatcherServletWrapper extends DispatcherServlet implements HttpSessionListener
{
	private static transient final Logger logger = LogManager.getLogger(DispatcherServletWrapper.class);
	private static final long serialVersionUID = 1L;

	public static final boolean profile = false;
	List<String> activeSessionIds = new ArrayList<String>();

	protected void doService(HttpServletRequest request, HttpServletResponse response) throws Exception 
	{
		response.setHeader("Content-type", "text/html");
		RequestStatistics requestStatistics = new RequestStatistics();
		requestStatistics.start();
		boolean noDB = request.getParameter("nodb") != null;
		TransportFactory.clearThreadTransport();
		ThreadScope.reset();
		ThreadScope.get().threadClassLoader = DispatcherServletWrapper.class.getClassLoader();
		Globals.clearMarshallingOptions();

		ThreadScope.get().fullDateFormat = new SimpleDateFormat("HH:mm, d MMM yy");
		ThreadScope.get().shortDateFormat = new SimpleDateFormat("d MMM yy");
		ThreadScope.get().requestURI = request.getRequestURI();

		if (!noDB)
		{
			Globals.startAllTransactions();
			if (profile)
			{
				Globals.session().createSQLQuery("set profiling_history_size = 100").executeUpdate();
				Globals.session().createSQLQuery("set profiling = 1").executeUpdate(); ////
			}
		}
		String uri = getRequestURI(request);
		if (uri.contains("longoperations/operationStatus") || uri.contains("getSessionData"))
			ThreadScope.get().disableLogging = true;

		if (OCHEMConfiguration.logLoad)
			LoadLogger.log();

		if(OCHEMConfiguration.verboseMode > 0){
			logger.info("[in] " + MemoryUtils.memorySummary());
			logger.info("[in] " + uri + " by " + getCurrentLogin(request, response));
		}

		try
		{
			super.doService(request, response);
			if (!noDB)
			{
				if (profile) // Midnighter on Jan 20, 2012
				{
					DecimalFormat df = new DecimalFormat("#.###");
					Globals.session().createSQLQuery("set profiling = 0").executeUpdate();
					List<Object[]> list = Globals.session().createSQLQuery("show profiles").list();

					for (Object[] row : list) {
						logger.info("[MySQL Profile] " + df.format(new Double("" + row[1])) + " sec:    " + row[2]);
					}
				}

				setStatus("Committing DB transactions...");
				Globals.commitAllTransactions();
			}
		}
		catch (UserFriendlyException e)
		{
			logger.error("Caught user friendly exception " + e.getMessage());
			Writer writer = response.getWriter();
			writer.write("Error: "+e.getMessage());
			writer.flush();
			writer.close();
			if (!noDB)
				Globals.rollbackAllTransactions();

		}
		catch (Exception e)
		{
			logger.error("Uncaught exception", e);
			// Midnighter
			// Thus, I believe, I resolve the problem,
			// described in
			// http://www.mikeschubert.com/archives/2006/08/javanetsocketex.html
			if (!noDB)
				try {
					Globals.rollbackAllTransactions();
				} catch (Exception e2) {
					logger.warn("Rollback while processing exception " + e.getMessage() + " failed.");
					// No probs. If rollback has failed, the original exception will be thrown in the code below.
				}
			throw e;
		}
		finally
		{
			setStatus("Finished");
			ThreadScope.get().threadFinished.fire();
			requestStatistics.stop();
			if(OCHEMConfiguration.verboseMode > 0){
				logger.info("[out] " + requestStatistics);
				logger.info("[out] " + MemoryUtils.memorySummary());
				logger.info("[out] " + uri + " by " + getCurrentLogin(request, response));
			}
		}
		//printPoolInfo();
	}

	private String getCurrentLogin(HttpServletRequest request, HttpServletResponse response)
	{
		if (response.isCommitted())
			return "";

		Session session = (Session) request.getSession().getAttribute("user-session");
		if (session != null)
			if (session.user == null)
				return "guest from " + ThreadScope.resolveRemoteAddr(request);
			else
				return session.user.login;
		else
			return "anonym";
	}

	public void init(ServletConfig config) throws ServletException 
	{
		super.init(config);
		//N.B. All logging is enabled
		
		//LogManager.getRootLogger().getAppender("console").addFilter(new ThreadLoggingFilter());
		try
		{
			String rootDir = config.getServletContext().getRealPath("/");
			Environment.rootDir = rootDir.substring(0, rootDir.length() - 1);
			Globals.initialize();
			Globals.jaxbContext = Globals.createJAXBContext();
		} catch (Exception e)
		{
			throw new ServletException(e);
		}
		logger.info("Dispatcher Servlet has been initialized with directory "+Globals.commonUploadDirectory);

	}

	//TODO -- change to proper Hibernate code

	private void cleaningJobs(){
		Globals.startAllTransactions();
		logger.info("Cleaning PendingTasks from tasks of guest and test users");
		Globals.session().createSQLQuery("delete PendingTask from PendingTask natural join Session "+
				"left join User using (user_id) where login is NULL or login = \"test\"").executeUpdate();
		Globals.commitAllTransactions();
		cleaningAnonymousModels();
		cleaningLostModels();
	}

	private void cleaningAnonymousModels(){
		Globals.startAllTransactions();
		logger.info("Cleaning Model from models of guest and test users");
		List <Integer>ids = Globals.session().createSQLQuery("select model_id from Model natural join Session " + 
				"left join User using (user_id) where login is NULL or login = \"test\"").list();
		Model.deleteModelsByIds(ids);
		Globals.commitAllTransactions();
	}

	private void cleaningLostModels(){
		Globals.startAllTransactions();
		logger.info("Cleaning Model from lost models");
		List <Integer>ids = Globals.session().createSQLQuery("select distinct Model.model_id from Model " + 
				"left join PendingTask using (task_id) where Model.task_id is not NULL and PendingTask.task_id is null").list();
		Model.deleteModelsByIds(ids);
		Globals.commitAllTransactions();
	}

	public void startupTest() throws ServletException
	{
		PeriodicTestRunner.loader = this.getClass().getClassLoader();	

		if (OCHEMConfiguration.testing)
		{
			PeriodicTestRunner runner = new PeriodicTestRunner();
			runner.filter.disableScheduler = true;
			new Thread(runner).start();
		}

		//This is test to inchi whether its working or not, if inchi is not working the we should stop the data base
		Globals.startAllTransactions();

		Criteria c = Globals.session().createCriteria(Molecule.class)
				.createCriteria("mapping2")
				.add(Restrictions.like("inchi2", "______________"))
				.setMaxResults(1);

		List<Molecule> molecules = c.list();
		if(molecules.size()!=0){
			Molecule molecule = molecules.get(0);
			logger.info("Got molecule with molid "+molecule.id+" and inchi "+molecule.mapping2.inchi2+" as or test victim");
			molecule.updateInchi();
			if(molecule.mapping1.inchi1.length() != 14)
				throw new ServletException("Database is down because Inchi is not working");
		}

		Globals.commitAllTransactions();
	}

	public void sessionCreated(HttpSessionEvent e) 
	{
		logger.info("Session has been created");
		activeSessionIds.add(e.getSession().getId());
		List<Session> session = Globals.session().createCriteria(Session.class).add(Restrictions.eq("session_string_id", e.getSession().getId())).list();
		if(session.size() > 0)
		{
			Session ses = session.get(0);
			ses.session_expired = 0;
		}
	}

	public void sessionDestroyed(HttpSessionEvent e) 
	{
		logger.info("Session has been destroyed");
		activeSessionIds.remove(e.getSession().getId());
		List<Session> session = Globals.session().createCriteria(Session.class).add(Restrictions.eq("session_string_id", e.getSession().getId())).list();
		if(session.size() > 0)
		{
			Session ses = session.get(0);
			ses.session_expired = 1;
		}
	}

	private String getRequestURI(HttpServletRequest request)
	{
		String user = "";
		if (request.getParameter("nodb") == null && Globals.userSession() != null && Globals.userSession().user != null)
			user = " by " + Globals.userSession().user;

		Enumeration<String> names = request.getParameterNames();
		List<String> namevaluepairs = new ArrayList<String>();
		while (names.hasMoreElements())
		{
			String name = names.nextElement();
			String[] values = request.getParameterValues(name);
			if (values == null)
				continue;
			if (values.length == 1)
				namevaluepairs.add(name+"="+values[0]);
			else
				namevaluepairs.add(name+"=["+StringUtils.join(values, ",")+"]");
		}
		return request.getRequestURI()+"?"+StringUtils.join(namevaluepairs, "&")+user;
	}

	private void setStatus(String status) {
		if (ThreadScope.get().operation != null)
			ThreadScope.get().operation.setStatus(status);
	}


}

/*
// A Log4J filter to disable logging for a particular thread
class ThreadLoggingFilter extends Filter
{
	@Override
	public int decide(LoggingEvent event)
	{
		return ThreadScope.get().disableLogging ? DENY : ACCEPT;
	}
}

*/
