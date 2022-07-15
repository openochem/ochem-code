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

import java.io.BufferedWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Session;
import org.hibernate.SessionFactory;
import org.hibernate.cfg.Configuration;

import qspr.metaserver.protocol.Command;

import com.eadmet.utils.config.GlobalConfigurator;
import com.eadmet.utils.mailer.Mailer;

// Abstract Hibernate-enabled server / Midnighter
// Can execute abstract commands
// Can report its status 

public abstract class AbstractServer 
{
	private static transient final Logger logger = LogManager.getLogger(AbstractServer.class);
	protected abstract Command executeCommand(Command request) throws Exception;
	protected abstract Configuration configureDbConnection(Configuration hibernateConf) throws Exception;
	protected abstract void serveGETRequest(HttpServletRequest req, HttpServletResponse res) throws Exception;

	private SessionFactory sessionFactory;
	protected Configuration dbConfiguration;
	protected ThreadLocal<HttpServletRequest> request = new ThreadLocal<HttpServletRequest>();
	public static AtomicInteger connections = new AtomicInteger(0);
	public static boolean debug = false;

	static {
		try {
			String[] configurations = {"/etc/ochem/metaserver.cfg"};
			GlobalConfigurator.configure(configurations);
			StringWriter sw = new StringWriter();
			GlobalConfigurator.export(new BufferedWriter(sw));
			logger.info("========= Actual final configuration ===========");
			logger.info(sw.toString());
			logger.info("========= Actual final configuration end========");
			Mailer.serverSignature = "Metaserver/Updateserver";
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	public static List<Object[]> allowedAddresses = new ArrayList<Object[]>();

	public static boolean checkAuth(String remoteAddress)
	{
		int i = 0;
		while (i < allowedAddresses.size()) 
		{
			Object[] addressTuple = AbstractServer.allowedAddresses.get(i);
			String address = (String)addressTuple[0];
			Date expiry = (Date)addressTuple[1];
			if (remoteAddress.startsWith(address))
				if (Calendar.getInstance().getTime().before(expiry))
					return true;
				else
				{
					allowedAddresses.remove(i);
					i--;
				};
				i++;
		}
		return false;
	}

	public AbstractServer(String id)
	{
		System.out.print("\nStarting server: " + id + "\n");
		try
		{
			logger.info("Trying to load " + id+".hibernate.cfg.xml");
			dbConfiguration = configureDbConnection(new Configuration().configure(id+".hibernate.cfg.xml"));
			sessionFactory = dbConfiguration.buildSessionFactory();
			logger.info("Database connection for "+id+" has been configured");
		}
		catch (Exception e)
		{
			logger.warn("Exception", e);
			logger.info("No database connection configured");
		}
	}

	public AbstractServer setRequest(HttpServletRequest request)
	{
		this.request.set(request);
		return this;
	}

	public static void debug(String st)
	{
		if (debug)
			logger.info(st);
	}



	public Session session()
	{	
		return sessionFactory.getCurrentSession();
	}
}
