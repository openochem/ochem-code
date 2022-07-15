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

package qspr.metaserver;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import javax.servlet.ServletConfig;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.protocol.Command;
import qspr.updateserver.UpdateServer;

import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.config.ConfigurableClass;
import com.eadmet.utils.config.ConfigurableProperty;
import com.eadmet.utils.mailer.Mailer;

@ConfigurableClass(name = "metaserver")
public class ServerServlet extends HttpServlet 
{
	private static transient final Logger logger = LogManager.getLogger(ServerServlet.class);

	@ConfigurableProperty(name = "allowed_ips")
	private static String allowedIPs = "127.0.0.1";

	public static final long serialVersionUID = 0x593F69EE;
	AbstractServer metaServer;
	AbstractServer updateServer;
	List<Integer> serverLoad = new ArrayList<Integer>();

	// For load statistics
	protected static AtomicInteger connectionAttempts = new AtomicInteger(0);
	protected static AtomicInteger acceptedConnectionAttempts = new AtomicInteger(0);
	private static int lastLoggedConnectionAttemts = 0;
	private static int lastLoggedAcceptedConnectionAttemts = 0;
	private long connectionAttemptsReset;

	public void doPost(HttpServletRequest req, HttpServletResponse res) throws ServletException, IOException
	{
		connectionAttempts.addAndGet(1);
		if (AbstractServer.connections.get() > 50)
		{
			logger.info("Too many connections, please retry later");
			res.sendError(HttpServletResponse.SC_SERVICE_UNAVAILABLE,"Too many connections, please retry later");
			return;
			//			throw new ServletException(new IOException("Too many connections, please retry later"));
		}
		Integer quenumber = AbstractServer.connections.addAndGet(1);
		acceptedConnectionAttempts.addAndGet(1);
		try
		{
			ObjectInputStream ois = new ObjectInputStream(req.getInputStream());
			if (AbstractServer.checkAuth(resolveRemoteAddr(req)))
			{
				Command command = (Command) ois.readObject();
				AbstractServer.debug("Entered with queue number execution <" + command + "> by "+command.senderId + " " + quenumber+" with id "+ Thread.currentThread().getId());
				Command responseCommand = getServer(req).setRequest(req).executeCommand(command);
				ObjectOutputStream oos = new ObjectOutputStream(res.getOutputStream());
				oos.writeObject(responseCommand);
				AbstractServer.debug("///Exit with queue number execution <" + command + "> by "+command.senderId + " " + quenumber+" with id "+ Thread.currentThread().getId());
			} else
			{
				logger.warn("Attempted access from an unallowed IP address " + resolveRemoteAddr(req));
				res.sendError(HttpServletResponse.SC_FORBIDDEN, "No access allowed from your machine " + resolveRemoteAddr(req));
			}

		} catch (Exception e)
		{
			logger.info("Here comes exception from IP address "+resolveRemoteAddr(req));
			logger.info(e);
			e.printStackTrace();
			res.sendError(HttpServletResponse.SC_INTERNAL_SERVER_ERROR,e.getMessage());
			//			throw new ServletException(e);
		}
		finally
		{
			AbstractServer.connections.addAndGet(-1);
		}
	}

	public void doGet(HttpServletRequest req, HttpServletResponse res) throws ServletException, IOException
	{
		try
		{
			if (req.getParameter("debug") != null)
				AbstractServer.debug = req.getParameter("debug").equals("true");

			if (AbstractServer.checkAuth(resolveRemoteAddr(req))) {
				logger.info("calling: "+ resolveRemoteAddr(req));
				AbstractServer server = getServer(req);
				logger.info("using: "+ server + " " + req + " " + res);
				server.serveGETRequest(req, res);
			}
			else
			{
				logger.warn("Attempted access from an unallowed IP address " + resolveRemoteAddr(req));
				res.sendError(HttpServletResponse.SC_FORBIDDEN, "No access allowed from your machine " + resolveRemoteAddr(req));
			}
		} catch (Exception e)
		{
			logger.info("Here comes exception from IP address "+ resolveRemoteAddr(req));
			logger.info(e);
			logger.info(e.getStackTrace());
			res.sendError(HttpServletResponse.SC_INTERNAL_SERVER_ERROR,e.getMessage());
		}
	}

	public void init(ServletConfig config)  throws ServletException
	{
		logger.info("Starting server servlet...");
		try
		{
			// Create actual servers
			metaServer = new MetaServer(MetaServer.METASERVER);
			updateServer = new UpdateServer("updateserver");
			String[] allowedAddresses = allowedIPs.split(",");
			logger.info("Allowed IP addresses: " + allowedIPs);
			AbstractServer.allowedAddresses.clear();
			for (String address : allowedAddresses)
			{
				Calendar expiry = Calendar.getInstance();
				expiry.add(Calendar.YEAR, 10);
				Object[] tuple = {address.trim(), expiry.getTime()};
				AbstractServer.allowedAddresses.add(tuple);
			}

			// Count average load per second
			new Thread(){

				public void run() {
					try
					{
						while (true)
						{
							Thread.sleep(5000);
							long time = Calendar.getInstance().getTimeInMillis()/1000;
							long total = ((connectionAttempts.get() - lastLoggedConnectionAttemts) / (time - connectionAttemptsReset));
							long accepted = ((acceptedConnectionAttempts.get() - lastLoggedAcceptedConnectionAttemts) / (time - connectionAttemptsReset));
							long lost = total - accepted;
							logger.info("Average connections per second: " + 
									total + " total, " + accepted + " accepted"+(lost > 0 ? ", " + lost + " lost" : ""));
							logger.info(MemoryUtils.memorySummary());
							checkMemory();
							lastLoggedConnectionAttemts = connectionAttempts.get();
							lastLoggedAcceptedConnectionAttemts = acceptedConnectionAttempts.get();
							connectionAttemptsReset = time;
						}
					}
					catch (Exception e)
					{
						logger.info(e);
					}
				}}.start();

				logger.info("Servlet server is running.");
		} catch (Exception e)
		{
			throw new ServletException(e);
		}
	}

	/**
	 * Notify developers on critical RAM shortage
	 */
	private boolean memoryDeficiencyMode = false;
	private void checkMemory()
	{
		double fraction = MemoryUtils.getCurrentMemoryUsedFraction();
		if (fraction > 0.80 && !memoryDeficiencyMode)
		{
			Mailer.notifyDevelopers("Metaserver memory usage more than 90%!", "Memory summary:\n" + MemoryUtils.memorySummary());
			memoryDeficiencyMode = true;
		}
		if (fraction < 0.8 && memoryDeficiencyMode)
		{
			Mailer.notifyDevelopers("Metaserver memory usage back to normal", "Memory summary:\n" + MemoryUtils.memorySummary());
			memoryDeficiencyMode = false;
		}
	}

	//	public boolean checkAuth(HttpServletRequest req, HttpServletResponse res) throws IOException
	//	{
	//		
	//		String remoteAddress = Globals.resolveRemoteAddr(req);
	//		for (String address : AbstractServer.allowedAddresses) 
	//		{
	//			if (remoteAddress.startsWith(address))
	//				return true;
	//		}
	//		res.sendError(HttpServletResponse.SC_FORBIDDEN, "No access allowed from your machine " + Globals.resolveRemoteAddr(req));
	//		return false;
	//	}

	public AbstractServer getServer(HttpServletRequest req)
	{
		if (req.getRequestURI().contains("update"))
			return updateServer;
		else
			return metaServer;
	}

	//Required to resolve IP when used under proxy by NGINX

	public static String resolveRemoteAddr(HttpServletRequest request) {
		String s = request.getHeader("X-Forwarded-For");
		if(s != null) return s;
		return request.getRemoteAddr();
	}

}
