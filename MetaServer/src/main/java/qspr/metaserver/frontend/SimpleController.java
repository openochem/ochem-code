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

package qspr.metaserver.frontend;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.servlet.ServletException;
import javax.servlet.ServletOutputStream;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.Transaction;
import org.hibernate.criterion.Order;

import qspr.metaserver.ArchivedTask;
import qspr.metaserver.MetaServer;
import qspr.metaserver.OnlinePeer;
import qspr.metaserver.StatisticsLog;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.NoSQLReference;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.NoSqlTransport;

import com.eadmet.utils.CryptUtils;


/**
 * Emulation of a simple and lightweight MVC framework for Metaserver
 * @author Midnighter
 *
 */

@SuppressWarnings("unchecked")
public class SimpleController 
{
	private static final Logger logger = LogManager.getLogger(SimpleController.class);
	public void invoke(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		logger.info("Invoking operation " + request.getParameter("action"));
		MetaServer ms = MetaServer.getInstance();
		synchronized (ms) 
		{
			Transaction tx = ms.session().beginTransaction();
			try
			{
				this.getClass().getMethod(request.getParameter("action"), HttpServletRequest.class, HttpServletResponse.class).invoke(this, request, response);
				tx.commit();
			}
			catch (Exception e)
			{
				e.printStackTrace();
				logger.warn(e);
				tx.rollback();
				throw e;
			}
		}
	}

	public void dev(HttpServletRequest request, HttpServletResponse response) throws IOException
	{
		response.setHeader("Content-type", "text/plain");
		PrintWriter writer = new PrintWriter(response.getOutputStream());
		writer.println(MetaServer.getInstance().onlinePeers.keySet());

		for (String  sId : MetaServer.getInstance().onlinePeers.keySet())
		{
			OnlinePeer peer = MetaServer.getInstance().onlinePeers.get(sId);
			writer.println("" + sId + " " + peer.getPing() + " " + peer.serverInfo);
		}
		writer.flush();
		writer.close();
	}

	public void getStatistics(HttpServletRequest request, HttpServletResponse response) throws IOException
	{
		int interval = 60;
		if (request.getParameter("interval") != null)
			interval = Integer.valueOf(request.getParameter("interval"));
		List<StatisticsLog> logs = MetaServer.getInstance().session().createCriteria(StatisticsLog.class)
				.addOrder(Order.desc("id"))
				.setMaxResults(interval)
				.list();

		response.setHeader("Content-type", "text/html");
		response.setContentType("text/html");
		PrintWriter writer = new PrintWriter(response.getOutputStream());
		writer.print("{\"logs\": [");
		for (StatisticsLog log : logs) {
			writer.println(String.format("{\"newTasks\": \"%s\", \"completedTasks\": \"%s\", \"assignedTasks\": \"%s\", \"connectionsPerSecond\": \"%s\", \"lostConnectionsPerSecond\": \"%s\", \"errors\": \"%s\", \"assignedTasks\": \"%s\"},", 
					log.newTasks, log.completedTasks, log.assignedTasks, log.connectionsPerSecond, log.lostConnectionsPerSecond, log.errors, log.assignedTasks));
		}

		writer.println("{}]}");
		writer.flush();
	}

	public void queue(HttpServletRequest request, HttpServletResponse response) throws IOException
	{
		response.setHeader("Content-type", "text/plain");
		response.setContentType("text/plain");

		PrintWriter writer = new PrintWriter(response.getOutputStream());
		try
		{
			if (MetaServer.getInstance().stopTaskAssignment)
				return;
			List<Object[]> rows = MetaServer.getInstance().session().createSQLQuery("select task_type, min_required_memory, count(*), max(priority) from Task where status='init' group by task_type, min_required_memory").list();
			for (Object[] objects : rows)
			{
				String taskType = objects[0].toString();
				Long minRequiredMemory = (objects[1] == null) ? 1 : Math.round(Double.valueOf(objects[1].toString()) / 1024); //In GB
				String taskCount = objects[2].toString();
				Long priority = (objects[1] == null) ? -100 : Math.round(Double.valueOf(objects[3].toString()));
				writer.println(taskType+","+minRequiredMemory+","+taskCount+","+priority);
			}
		}
		finally
		{
			writer.flush();
		}
	}

	private void getFile(String fileName, HttpServletResponse response) throws IOException, Exception
	{
		if (!fileName.matches("[A-Za-z\\-]+\\.[A-Za-z]+"))
			throw new ServletException("Not allowed file name: " + fileName);

		StringWriter sWriter = new StringWriter();
		PrintWriter writer = new PrintWriter(sWriter);
		File versionTemplate = new File("/etc/ochem/" + fileName);
		if (versionTemplate.exists())
		{
			BufferedReader reader = new BufferedReader(new FileReader(versionTemplate));
			String st;
			while ((st = reader.readLine()) != null)
				writer.println(st);
			reader.close();

			OutputStream os = response.getOutputStream();
			os.write(CryptUtils.desEncode(sWriter.toString()).getBytes());
			os.close();
		}
		else
			throw new ServletException("The configuration file " + fileName + " could not be found");
		writer.close();
	}

	public void versionTemplate(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		getFile(Command.VERSIONTEMPLATEXML, response);
	}

	public void getFile(HttpServletRequest request, HttpServletResponse response) throws Exception
	{
		getFile(request.getParameter("file"), response);
	}

	public void charts(HttpServletRequest req, HttpServletResponse res) throws ServletException, IOException
	{
		req.getRequestDispatcher("/jsp/charts.jsp").forward(req, res);
	}


	private void debugfilesList(HttpServletRequest req, HttpServletResponse res) throws IOException, ServletException
	{
		List<Map<String, String>> files = NoSqlTransport.listFiles("debug", "files");
		Map<Integer, List<Map<String, String>>> resultt = new HashMap<Integer, List<Map<String, String>>>();
		Map<String, List<Map<String, String>>> results = new HashMap<String, List<Map<String, String>>>();
		for (Map<String, String> file : files) 
		{
			if (file.get("task_id") != null)
			{
				Integer task_id = Integer.valueOf(file.get("task_id"));
				List<Map<String, String>> taskList = resultt.get(task_id);
				if (taskList == null)
					resultt.put(task_id, taskList = new ArrayList<Map<String, String>>());
				taskList.add(file);
			} else
				if (file.get("server_name") != null)
				{
					String server_name = file.get("server_name");
					List<Map<String, String>> taskList = results.get(server_name);
					if (taskList == null)
						results.put(server_name, taskList = new ArrayList<Map<String, String>>());
					taskList.add(file);				
				}
		}
		req.setAttribute("files_tasks", resultt);
		req.setAttribute("files_servers", results);
		req.getRequestDispatcher("/jsp/debugfiles.jsp").forward(req, res);
	}

	private void debugfilesDownload(HttpServletRequest req, HttpServletResponse res) throws ServletException, IOException
	{
		String id = req.getParameter("id");
		Map<String, String> file = NoSqlTransport.listFile("debug", "files", id);
		if (file == null)
			return;
		byte[] data = NoSqlTransport.getDataSafely(new NoSQLReference(id, "debug", "files"));
		res.setHeader("Content-Disposition", "attachment; filename="+file.get("filename"));
		ServletOutputStream sos = res.getOutputStream();
		sos.write(data);
		sos.flush();
		sos.close();
	}

	private void debugfilesDelete(HttpServletRequest req, HttpServletResponse res) throws ServletException, IOException
	{
		String id = req.getParameter("id");
		Map<String, String> file = NoSqlTransport.listFile("debug", "files", id);
		if (file == null)
			return;
		NoSqlTransport.deleteFile("debug", "files", id);

		res.sendRedirect("?action=debugfiles");
	}

	private void debugfilesClear(HttpServletRequest req, HttpServletResponse res) throws ServletException, IOException
	{
		List<Map<String, String>> files = NoSqlTransport.listFiles("debug", "files");
		for (Map<String, String> file : files) 
		{
			String id = file.get("_id");
			NoSqlTransport.deleteFile("debug", "files", id);
		}
		res.sendRedirect("?action=debugfiles");
	}

	public void debugfiles(HttpServletRequest req, HttpServletResponse res) throws ServletException, IOException
	{
		String saction = req.getParameter("saction");
		if (saction == null || "list".equals(saction))
			debugfilesList(req, res);
		else if ("download".equals(saction))
			debugfilesDownload(req, res);
		else if ("delete".equals(saction))
			debugfilesDelete(req, res);
		else if ("clear".equals(saction))
			debugfilesClear(req, res);
	}	

	public void resettasks(HttpServletRequest req, HttpServletResponse res) throws ServletException, IOException
	{
		int minmem = 3000;	
		MetaServer.getInstance().session().createSQLQuery("update Task set min_required_memory="+minmem+", resubmitted=1 where resubmitted > 0  and min_required_memory > "+minmem+"").list();
		MetaServer.getInstance().taskQueue.clear();
		MetaServer.getInstance().session().createSQLQuery("update Task set min_required_memory="+minmem+", resubmitted=0 where resubmitted > 0  and min_required_memory > "+minmem+"").list();
		res.sendRedirect(MetaServer.getLocalURL(req));
	}	

	public void taskProfile(HttpServletRequest req, HttpServletResponse res)
	{
		try
		{
			Integer id = Integer.valueOf(req.getParameter("id"));
			Task t = (Task) MetaServer.getInstance().session().get(Task.class, id);

			TaskTreeNode root = null;
			if (t == null)
			{
				ArchivedTask at = (ArchivedTask) MetaServer.getInstance().session().get(ArchivedTask.class, id);
				root = new TaskTreeNode(at);
			}
			else
				root = new TaskTreeNode(t);

			req.setAttribute("rootNode", root);
			req.getRequestDispatcher("/jsp/taskprofile.jsp").forward(req, res);
		}
		catch (Exception e)
		{
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

	public void serversAvailability(HttpServletRequest req, HttpServletResponse res) throws ServletException, IOException
	{
		ServersAvailbilityData data = new ServersAvailbilityData();
		data.addPeers(MetaServer.getInstance().onlinePeers.values());

		List<Object[]> rows = MetaServer.getInstance().session().createSQLQuery("select task_type, priority from Task where status='init'").list();
		for (Object[] row : rows)
			data.addTask((String) row[0], (int) Math.round(Math.floor(Double.valueOf("" + row[1]))), false);

		for (Task task : MetaServer.getInstance().assignedTasks.values())
			data.addTask(task.taskType, (int) Math.round(Math.floor(task.priority)), true);

		req.setAttribute("data", data);
		req.getRequestDispatcher("/jsp/servers-availability.jsp").forward(req, res);
	}

	public void addTaskBinding(HttpServletRequest req, HttpServletResponse res) throws ServletException, IOException
	{
		MetaServer.getInstance().forcedTaskBindings.put(req.getParameter("taskType"), req.getParameter("server"));
	}

	public void removeTaskBinding(HttpServletRequest req, HttpServletResponse res) throws ServletException, IOException
	{
		MetaServer.getInstance().forcedTaskBindings.remove(req.getParameter("taskType"));
	}

	private static SimpleController instance;

	/**
	 * Singleton implementation
	 * @return The instance of the SimpleController
	 */
	public static SimpleController getInstance()
	{
		if (instance == null)
			instance = new SimpleController();
		return instance;
	}

	public static void main(String[] args)
	{
		System.out.println("AA.aa".matches("[A-Za-z\\-]+\\.[A-Za-z]+"));
	}
}
