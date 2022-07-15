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

package qspr;

import java.text.SimpleDateFormat;
import java.util.List;

import javax.servlet.http.HttpServletRequest;

import org.apache.commons.lang3.mutable.Mutable;
import org.apache.commons.lang3.mutable.MutableObject;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Transaction;
import org.springframework.web.multipart.MultipartHttpServletRequest;

import qspr.entities.ExperimentalProperty;
import qspr.metaserver.Event;
import qspr.util.Operation;

/**
 * Stores the variables related to the current thread context.
 * @author midnighter
 *
 */
@SuppressWarnings("rawtypes")
public class ThreadScope 
{
	public HttpServletRequest localRequest;
	public MultipartHttpServletRequest localMpRequest;
	public Mutable<Transaction> transaction = new MutableObject<Transaction>();
	public Mutable<Transaction> alternateTransaction = new MutableObject<Transaction>();
	public ClassLoader threadClassLoader;
	public qspr.entities.Session userSession;
	public Boolean exclusiveLockOnSession;
	public Boolean updateSessionTime;
	public Operation operation;
	public SimpleDateFormat fullDateFormat;
	public SimpleDateFormat shortDateFormat;
	public String controller;
	public String requestURI;
	public Boolean disableTrackChanges;
	public Boolean considerPredicatesInStatistics;
	public List<MarshallingOption> marshallingOptions;

	/**
	 * Flag to indicate disabled Log4J logging for this thread
	 */
	public boolean disableLogging = false;

	public String context; // Sometimes its good to know what is this thread about globally

	public Event threadFinished = new Event(this);
	public Event transactionCommit = new Event(this);
	public Event transactionRollback = new Event(this);

	public Event<ExperimentalProperty> recordApproved = new Event<ExperimentalProperty>(this);


	public static ThreadScope get()
	{
		ThreadScope scope = threadScope.get();
		if (scope == null)
			threadScope.set(scope = new ThreadScope());

		return scope;
	}

	public HttpServletRequest getHttpServletRequest()
	{
		if (localMpRequest != null)
			return localMpRequest;
		else
			return localRequest;
	}

	public static void reset()
	{
		threadScope.remove();
	}

	public static void setStatus(String status)
	{
		setStatus(status, logger);
	}

	public static void setStatus(String status, Logger logger)
	{
		if(OCHEMConfiguration.verboseMode > 0)
			logger.info(status);
		if (get().operation != null)
			get().operation.setStatus(status);
	}

	private static ThreadLocal<ThreadScope> threadScope = new ThreadLocal<ThreadScope>();

	private static final Logger logger = LogManager.getLogger(ThreadScope.class);

	//Required to resolve IP when used under proxy by NGINX

	public static String resolveRemoteAddr(HttpServletRequest request) {
		String s = request.getHeader("X-Forwarded-For");
		if(s != null) return s;
		return request.getRemoteAddr();
	}

}
