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

package qspr.util;

import java.util.HashSet;
import java.util.Set;

import javax.servlet.http.HttpServletRequest;

import org.apache.commons.lang3.mutable.Mutable;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Transaction;

import qspr.Environment;
import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.entities.PendingTask;
import qspr.entities.Session;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.Transport;
import qspr.metaserver.transport.TransportFactory;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.mailer.Mailer;

/**
 * A wrapper class that takes care of miscelleneous required initialisations (database connection, JAXB context, etc)  
 * @author midnighter
 *
 */
abstract public class WrapperThread extends Thread
{
	private static transient final Logger logger = LogManager.getLogger(WrapperThread.class);

	public Session userSession = Globals.userSession();
	public HttpServletRequest httpRequest = null;
	public Transport parentTransport = TransportFactory.getThreadTransport();
	public StatusTracker statusTracker = new StatusTracker(logger);

	/**
	 * The pending task associated with this thread (optional)
	 */
	public PendingTask pTask;


	public boolean updateSessionTime = true;
	public boolean cancelRequested = false;

	/**
	 * The operation associated with this thread (optional)
	 */
	public Operation operation;

	private ClassLoader parentClassLoader;
	public Exception exception;

	/**
	 * The list of all currently running WrapperThread ancestors
	 */
	public static Set<WrapperThread> runningThreads = new HashSet<WrapperThread>();

	public void run()
	{
		try
		{
			if (Globals.jaxbContext == null)
				Globals.jaxbContext = Globals.createJAXBContext();
		}
		catch (Exception e)
		{
			throw new RuntimeException(e);
		}

		setStatus("Initializing...");
		ThreadScope.get().updateSessionTime = updateSessionTime;
		ThreadScope.get().operation = operation;
		TransportFactory.setThreadTransport(parentTransport);
		runningThreads.add(this);
		try 
		{
			if (Environment.rootDir == null)
				Environment.rootDir = "/ews/Toxicity";
			if (ThreadScope.get().threadClassLoader == null)
				ThreadScope.get().threadClassLoader = parentClassLoader;
			ThreadScope.get().localRequest = httpRequest;

			if (!OCHEMConfiguration.disableDB)
			{
				Globals.startAllTransactions();
				if (userSession != null)
				{
					Session oldSession = userSession;
					ThreadScope.get().userSession = userSession = (Session)Globals.session().get(Session.class, userSession.id);
					if (userSession != null)
						userSession.copySessionInfoFrom(oldSession);
				}
			}
			else
				ThreadScope.get().userSession = userSession;
			wrapped();

			safeCommit(ThreadScope.get().transaction, "main(qspr)");
			safeCommit(ThreadScope.get().alternateTransaction, "alternate(fragments)");
			if (cancelRequested)
				setStatus("Finished - operation has been cancelled");
		} catch (Exception e) 
		{
			exception = e;
			setStatus("Error: "+e.getMessage());
			logger.info("[WrapperThread] ===== Catched Exception occured in WrapperThread =====" + this.getClass().getSimpleName());
			if (pTask != null)
				logger.info("Error happened when processing pending task " + pTask.id + ", metaserver task ID " + pTask.taskId);
			e.printStackTrace();
			Globals.rollbackAllTransactions();
		}
		finally
		{
			// The thread has finished - store its last status
			//if (operationId != null)
			//	operationStatuses.put(operationId, status);
			ThreadScope.reset();
			TransportFactory.clearThreadTransport();
			runningThreads.remove(this);
			System.gc();
		}
	}

	public void setStatus(String newStatus)
	{
		statusTracker.set(newStatus);
		if (operation != null)
			operation.setStatus(newStatus);
		else
			logger.info("operation is null: " + newStatus);
	}

	public String getStatus()
	{
		return statusTracker.get();
	}

	public WrapperThread()
	{
		try
		{
			parentClassLoader = ThreadScope.get().threadClassLoader;
			if (ThreadScope.get().localRequest != null && ThreadScope.get().localRequest.getParameter("operation-id") != null)
			{
				// Register an operation for this thread
				setOperationID(Long.valueOf(ThreadScope.get().localRequest.getParameter("operation-id")));
			}
			else
				System.out.println( ">>>> ThreadScope.get().localRequest " + ThreadScope.get().localRequest + " ThreadScope.get().localRequest.getParameter(\"operation-id\") " + ThreadScope.get().localRequest.getParameter("operation-id"));
		}
		catch (Exception e)
		{
			// When invoked form a web-service, "getParameter" fails for an unknown reason. Ignore it: operations API does not apply to web-services.
		}
	}

	public WrapperThread setOperationID(long operationID)
	{
		setName("Operation-" + operationID);
		operation = new Operation(operationID);
		System.out.println(">>>> setOperationID:  starting new operation with id: " + operationID);
		return this;
	}

	public WrapperThread newOperation()
	{
		operation = new Operation(Operation.generateID());
		System.out.println(">>>> newOperation:  creating new operation with id: " + operation.operationId);
		return this;
	}

	public static WrapperThread getByOperationID(long operationID)
	{
		for (WrapperThread t : WrapperThread.runningThreads)
			if (t.operation != null && t.operation.operationId == operationID)
				return t;
		return null;
	}


	/**
	 * Commit without throwing exceptions.
	 * In case of a failure, log it and send email to developers 
	 * @param tx
	 * @param name
	 */
	private void safeCommit(Mutable<Transaction> transaction, String name)
	{
		Transaction tx = transaction.getValue();
		try
		{
			if (tx != null && tx.isActive())
				tx.commit();
		}
		catch (Exception e)
		{
			logger.error("Failed to commit a transaction " + name + "! may result into loss of data!");
			e.printStackTrace();

			this.exception = e;

			// This is quite bad. Notify developers so that they know about the problem. // Midnighter on Jan 24, 2012
			String msg = "We could not commit a transaction in " + this.getClass().getName() + "\nThis is bad and can result into inconsistencies.\n";
			if (pTask != null)
			{
				msg += "\nPending task " + pTask.id;
				if (pTask.model != null)
					msg += "\nModel: " + pTask.model.name;
			}

			msg += "\n\nStacktrace: \n" + OCHEMUtils.exceptionToString(e);

			Mailer.notifyDevelopers("Could not commit a transaction in WrapperThread", msg);
		} finally
		{
			transaction.setValue(null);
		}
	}

	public void registerPendingTask(PendingTask pTask)
	{
		pTask.status = Task.INIT;
		pTask.updatePostedTime();
		Globals.session().saveOrUpdate(pTask);
		this.pTask = pTask;
	}

	/**
	 * Has this thread terminated because of a failure?
	 * @return
	 */
	public boolean isFailed()
	{
		return exception != null;
	}

	/**
	 * Get a thread associated wihh a particular pending task
	 * @param pTaskId
	 * @return
	 */
	public static WrapperThread getPendingTaskThread(long pTaskId)
	{
		for (WrapperThread thread : runningThreads)
			if (thread.pTask != null && thread.pTask.id != null && thread.pTask.id.equals(pTaskId))
				return thread;

		return null;
	}


	/**
	 * Is this thread currently running?
	 * @return
	 */
	public boolean isRunning()
	{
		return runningThreads.contains(this);
	}

	public void start()
	{
		if (isRunning())
			throw new UserFriendlyException("Can't start a thread twice: " + this);
		runningThreads.add(this);
		super.start();
	}

	/**
	 * The actual functionality of the thread wrapped in a helper code (init, exception handling, etc)
	 * @throws Exception
	 */
	abstract public void wrapped() throws Exception;

}
