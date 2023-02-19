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

package qspr.metaserver.serv;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.lang.reflect.Method;
import java.net.InetAddress;
import java.net.MalformedURLException;
import java.net.URI;
import java.net.UnknownHostException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.core.config.Configurator;
import org.junit.runner.Description;
import org.junit.runner.JUnitCore;
import org.junit.runner.Request;
import org.junit.runner.Result;
import org.junit.runner.manipulation.Filter;
import org.junit.runner.notification.Failure;
import org.junit.runner.notification.RunListener;

import qspr.configuration.JAXBContextFactory;
import qspr.metaserver.CalculationServer;
import qspr.metaserver.ServerPool;
import qspr.metaserver.cs.util.CpuMeasurementTool;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.ServerInfo;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.DebugLevel;
import qspr.metaserver.serv.ServerRunnerConfiguration.LocalAppConfig;
import qspr.metaserver.tests.DescriptorServerTest;
import qspr.metaserver.tests.MainTest;
import qspr.metaserver.tests.TaskTest;
import qspr.metaserver.transport.CSTransport;
import qspr.metaserver.transport.DataReferenceFactory;
import qspr.metaserver.transport.NoSqlTransport;
import qspr.workflow.utils.QSPRConstants;
import qspr.tests.ReportGeneratorTestListener;
import qspr.workflow.utils.ProcessUtils;

import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.exe.OutReader;

/**
	 MultiServer runs multiple types of tasks, one task at a time
	 Runs tests on startup.
	 Requires a valid server home directory.

	 @author midnighter
 */

public class MultiServer {
	// Declarations of constants
	static final long MINIMUM_DISK_SPACE = 500 * 1024 * 1024; // 500 MB
	static final long CRITICAL_DISK_SPACE = 1000 ; // ca 1GB in MB
	static final long TEST_PASSED_RERUN_TIME = 240 * 60 * 60 * 1000; // 10 days
	static final long TEST_FAILED_RERUN_TIME = 60 * 60 * 1000; // 1 hour
	static final long MAXIMUM_TIME_METASERVER_UNAVAILABLE = 12 * 60 * 60 * 1000; // 12 hours

	static final long testExpieryRandom = Math.round(Math.random() * TEST_FAILED_RERUN_TIME); // not to have all tests running simultaneously
	private long testsExpieryTime = TEST_FAILED_RERUN_TIME; // by default, time to rerun for failed tasks

	private String workingDirectory = ".";
	private CSTransport transport = new CSTransport();
	private boolean runIt;
	private String sid = null;
	private static String version = null;
	private CalculationServer currentCalculationServer;
	private ServerInfo serverInfo = new ServerInfo();
	private boolean finishAfterDone = true;
	private boolean generateTestReport = false;
	// Tests control
	private long lastTestsRun;
	private boolean restartRequested = false;
	private boolean terminateRequested = false;
	private boolean clean = true;
	private String configurationXml = null;
	public boolean skipTests = false;

	private static TunnelPrintWriter out = null;
	private static CalculationThread calculationThread = null;

	private JAXBContext jaxbContext;

	CpuMeasurementTool cpuMeasurementTool = new CpuMeasurementTool();

	private long idleSince; // Since when the server has been running idle

	public static void main(String[] args) throws InterruptedException, IOException {
		try {
			MultiServer m = new MultiServer();
			ServerRunnerConfiguration c = m.initServer(args);
			ServerRunnerConfiguration.instance = c;
			m.run(c);
		} catch (MalformedURLException e) {
			notifyRunner(ExchangeCommands.MULTISERVER_TERMINATE, "Metaserver URL is invalid");
		} catch (Exception e) {
			if (out != null)
				e.printStackTrace(out);
			notifyRunner(ExchangeCommands.MULTISERVER_TERMINATE, "Exception when running server: " + e.getMessage());
		} catch (Throwable e) {
			if (out != null)
				e.printStackTrace(out);
			notifyRunner(ExchangeCommands.MULTISERVER_TERMINATE, "Critical Exception when running server: " + e.getMessage());
		}

	}


	public void run(ServerRunnerConfiguration configuration) throws MalformedURLException, ClassNotFoundException, InterruptedException {
		// Main working loop of the server
		runIt = true;
		idleSince = Calendar.getInstance().getTimeInMillis();

		boolean justStarted = true;

		long lastGC = Calendar.getInstance().getTimeInMillis() + 3*60*1000; // 3 minutes delay to avoid problem with overload for running of tests

		serverInfo.usedMemory =   MemoryUtils.getCurrentMemoryUsage();
		serverInfo.peakMemory = MemoryUtils.getPeakMemoryUsage();


		while (runIt)
			try {

				if(Calendar.getInstance().getTimeInMillis() - lastGC > 60*1000){
					System.gc(); // our goal is allocate as minimum memory as possible and we do it once per minute
					lastGC = Calendar.getInstance().getTimeInMillis(); 
					serverInfo.usedMemory =   MemoryUtils.getCurrentMemoryUsage();
					serverInfo.peakMemory = MemoryUtils.getPeakMemoryUsage();
				}

				if (calculationThread == null) {
					if (idleMode(configuration)) {
						// if we are still in idle mode, we will sleep...
						if (finishAfterDone) // corresponds to traditional mode -- server finishes after calculation
							notifyRunner(ExchangeCommands.MULTISERVER_PING, null); // modify server runner that we are in idle mode and are sleeping ...
						if (!justStarted)
							Thread.sleep(configuration.sleepTime * 1000); // avoid the delay after the starting server
						justStarted = false;
						if (configuration.maximumIdleTime != null
								&& (Calendar.getInstance().getTimeInMillis() - idleSince > configuration.maximumIdleTime * 1000)) {
							out.println("Exceeded the maximum idle time. The server is scheduled for terminating.");
							terminateRequested = true;
						}
					}
				} else if (!calculationThread.isAlive()) {
					onTaskFinished(null);
					idleSince = Calendar.getInstance().getTimeInMillis();
				} else {
					taskCalculationMode();
					synchronized (calculationThread) {
						if (calculationThread.isAlive())
							calculationThread.wait(getSleepTimeSeconds(calculationThread.getRunningTimeSeconds(), configuration.sleepTimeWhileCalculating) * 1000);
					}
				}
			} catch (IOException e) {
				serverInfo.configurationXml = configurationXml; // Since we may need to re-register us with the server
				long unavailability = (Calendar.getInstance().getTimeInMillis() - transport.lastSuccessfullConnect);
				out.println("Metaserver is unavailable for " + unavailability / 1000 + " seconds: " + e.getMessage());
				if (unavailability > MAXIMUM_TIME_METASERVER_UNAVAILABLE) // more than 12 hours
				{
					runIt = false;
					out.println("Terminating, because metaserver is not available for more than 12 hours");
				}
			}
	}

	/**
	 * Progressive sleep time.
	 * The idea: longer tasks do not need to be as interactive as short tasks
	 */
	private long getSleepTimeSeconds(long runningTimeSeconds, long maxSleepTime) {
		long sleep = maxSleepTime;
		if (runningTimeSeconds < 60)
			sleep = 5;
		else if (runningTimeSeconds < 300)
			sleep = 10;
		else if (runningTimeSeconds < 20 * 60)
			sleep = 15;
		else
			return maxSleepTime;

		return Math.min(sleep, maxSleepTime);
	}

	public void executeUserCommand(String command) {
		command = command.toUpperCase();
		if (command.equals("Q"))
			runIt = false;
		else if (command.equals(ExchangeCommands.SERVER_RESTART_COMMAND))
			restartRequested = true;
		else if (command.equals(ExchangeCommands.SERVER_TERMINATE_COMMAND))
			terminateRequested = true;
		else if (command.startsWith("SETPARAM"))
		{
			// A command of the format: SETPARAM <APPNAME> <PARAMNAME> <VALUE>
			String[] parts  = command.split(" ");
			ServerPool.getInstance().getServer(parts[1], false).setParam(parts[2], parts[3]);
			ServerRunnerConfiguration.instance.setApplicationParam(parts[1], parts[2], parts[3]);

			saveConfiguration(ServerRunnerConfiguration.instance);
		}

		out.println("Executing command <" + command + ">");
	}

	private void saveConfiguration(ServerRunnerConfiguration config)
	{
		try
		{
			Marshaller marshaller = jaxbContext.createMarshaller();
			marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, new Boolean(true));
			marshaller.marshal(config, new File(workingDirectory + "/" + ServerRunner.VERSIONXML));
		} catch (JAXBException e)
		{
			e.printStackTrace();
		}
	}

	/**
	 * Provides checking of initial status before running calculations; creates directory and configuration
	 * @throws Exception 
	 */
	private ServerRunnerConfiguration initServer(String[] args) throws Exception {


		if (args.length == 0)
			notifyRunner(ExchangeCommands.MULTISERVER_TERMINATE, "Server home directory not specified!");

		for (String arg : args)
			if (arg.equals("--skip-tests"))
				skipTests = true;

		transport.deepSleepTime = 10 * 1000; // We use a sparing deep sleep interval here

		workingDirectory = args[0];
		out = new TunnelPrintWriter(new FileOutputStream(workingDirectory + "/output/server.out", true), null, true, false);
		redirectLog4jTo(out);
		if (args.length > 1 && args[1].equals("--standalone")) {
			finishAfterDone = false;
			new File(workingDirectory + ExchangeCommands.MULTISERVER_TESTS_FILE_NAME).delete();
		}

		if (!isEnoughDiskSpace(MINIMUM_DISK_SPACE))
			notifyRunner(ExchangeCommands.MULTISERVER_TERMINATE, "There is not enough disk space.");

		jaxbContext = JAXBContextFactory.get("qspr.metaserver.serv");
		Unmarshaller unmarshaller = jaxbContext.createUnmarshaller();
		ServerRunnerConfiguration configuration = CalculationServerConfigurator.readConfigurationFile(jaxbContext, workingDirectory);

		generateTestReport = (configuration.generateTestReport == null) ? false : configuration.generateTestReport;

		if (configuration.metaserverURL != null && !configuration.metaserverURL.equals("")) {
			configuration.metaserverURL = CSTransport.addSeparator(configuration.metaserverURL);
			transport.serverURL = CSTransport.defaultServerURL = configuration.metaserverURL;
		}
		if (configuration.mongoDbURL != null && !configuration.mongoDbURL.equals(""))
			NoSqlTransport.host = configuration.mongoDbURL;
		out.println("MongoDB URL: " + NoSqlTransport.host);
		out.println("Metaserver URL: " + CSTransport.defaultServerURL);
		out.println("OCHEM URL (for web-service invocations): " + configuration.ochemURL);

		int currentId = ProcessUtils.getCurrentPid();
		int currentOwnerId = ProcessUtils.getParentPid();

		// Check if we have a concurrent server running
		if (configuration.processId != null && configuration.processId != currentId
				&& ProcessUtils.isProcessRunning(configuration.processId, configuration.parentProcessId))
			notifyRunner(ExchangeCommands.MULTISERVER_TERMINATE, "There is another server already running from the same directory. Terminating...");

		// Save current process id to version.xml
		configuration.processId = currentId;
		configuration.parentProcessId = currentOwnerId;

		Applications applications = (Applications) unmarshaller.unmarshal(this.getClass().getClassLoader().getResourceAsStream("applications.xml"));
		ServerRunnerConfiguration.allApplications = applications;

		if (configuration.applications == null || configuration.applications.size() == 0)
			notifyRunner(ExchangeCommands.MULTISERVER_TERMINATE, "To complete setup you must provide correct " + ServerRunner.VERSIONXML +"  file");

		try {
			// FIX to work inside of the docker
			try{ // if "ochem-mongo" is available, we use it and also update path to mariadb
				InetAddress inet[] = InetAddress.getAllByName("ochem-mongo");
				if(inet!=null && inet.length>0) {
					configuration.mongoDbURL = NoSqlTransport.host = "mongodb://ochem-mongo";
					if(configuration.applicationConfigurations != null)
						for(LocalAppConfig ap: configuration.applicationConfigurations)
							for(ApplicationParam param: ap.params) {
								if(param.name.equals("DBURL")) {
									String cleanURI = param.value.substring(5);
									URI uri = URI.create(cleanURI);
									String host = uri.getHost();
									String rest=cleanURI.substring(cleanURI.indexOf(host)+host.length());
									param.value =  "jdbc:mariadb://ochem-mariadb"+rest;
								}
							}
				}
			}catch (Exception e) {}

			// trying whether MongoDB is available
			int trials = NoSqlTransport.trials;
			NoSqlTransport.trials = 0;
			DataReferenceFactory.createReferencer().saveReference("test", QSPRConstants.WEBSERVICE_DATABASE);
			NoSqlTransport.trials = trials;
		}catch(Exception e) {
			out.println("Mongodb is not available at: " + configuration.mongoDbURL);
			System.exit(1);
		}

		saveConfiguration(configuration);

		version = serverInfo.version = configuration.currentVersion;
		serverInfo.ipAddress = InetAddress.getLocalHost().getHostAddress();
		serverInfo.workingDirectory = workingDirectory;
		serverInfo.name = configuration.sid;
		serverInfo.minimumPriority = configuration.getMinimumPriority();
		serverInfo.availableMemory = configuration.memoryLimit;
		if (serverInfo.minimumPriority != null)
			out.println("Minimum priority: " + serverInfo.minimumPriority);
		serverInfo.random = Math.round(Math.random() * 10000);
		serverInfo.owner = ProcessUtils.getComputerOwner();
		configurationXml = serverInfo.configurationXml = getFileAsString(workingDirectory + "/" + ServerRunner.VERSIONXML);

		// Unique SID is combined from user-specified SID and working directory
		sid = MultiServer.getAutomaticSID(workingDirectory);
		ServerPool.getInstance().sid = sid;

		configureSupportedServers(applications, configuration);

		out.println("Multiserver started, session id is " + serverInfo.random + ". Supported tasks: " + serverInfo.supportedTaskTypes);
		new OutReader(System.in, Thread.currentThread()) {
			public void onLineRead(String line) {
				executeUserCommand(line);
			}
		}.start();

		return configuration;
	}

	/**
	 * There is request to restart servers
	 * 
	 * @return
	 */

	private boolean shoudRestart() {
		return finishAfterDone && !clean;
	}

	/**
	 * returns true if there was no message from the metaserver, i.e. we are still in idle mode
	 */

	public boolean idleMode(ServerRunnerConfiguration configuration) throws ClassNotFoundException, MalformedURLException, IOException, InterruptedException {
		serverInfo.status = "idle";

		// Constantly update the minimum priority, for it can be changed with time // Midnighter on Nov 16, 2011
		serverInfo.minimumPriority = configuration.getMinimumPriority();

		if (restartRequested || shoudRestart() || (Calendar.getInstance().getTimeInMillis() - lastTestsRun >= testsExpieryTime + testExpieryRandom && !skipTests)) {
			notifyRunner(ExchangeCommands.MULTISERVER_RESTART, "Restarting to rerun tests");
			return true;
		}

		if (!isEnoughDiskSpace(MINIMUM_DISK_SPACE)) {
			notifyRunner(ExchangeCommands.MULTISERVER_TERMINATE, "Terminate due to critical disk space.");
			return true;
		}

		if (terminateRequested) {
			notifyRunner(ExchangeCommands.MULTISERVER_TERMINATE, "Terminate requested.");
			return true;
		}
		serverInfo.setSupportedTasks(configuration.getSupportedTaskTypes(), configuration.currentVersion);
		serverInfo.cpuUsage = cpuMeasurementTool.getLoads()[2];
		Command response = transport.executeCommand(new Command(Command.CS_REGISTER, serverInfo).sid(sid));
		serverInfo.configurationXml = null; // we successfully reported our configuration

		if (response == null)
			return true;

		switch (response.id) {
		case Command.MS_TERMINATE:
			// If the metaserver requested us to terminate
			runIt = false;
			notifyRunner(ExchangeCommands.MULTISERVER_TERMINATE, "Terminating server upon request from the Metaserver.");
			return true;
		case Command.MS_RESTART:
			notifyRunner(ExchangeCommands.MULTISERVER_RETEST, "Restarting upon request from the Metaserver.");
			restartRequested = true;
			return true;
		case Command.MS_COMMAND:
			executeUserCommand(response.data.toString());
			return true;
		case Command.MS_UPDATEREQUIRED:
			notifyRunner(ExchangeCommands.MULTISERVER_UPDATEREQUIRED, "Notifying the upsteram runner that an update is required");
			return true;
		case Command.MS_GET_LOGS:
			currentCalculationServer.takeLogsSnapshot(null);
			sendDebugFiles(null);
			return true;
		default:
			// We have a task to calculate!
			notifyRunner(ExchangeCommands.MULTISERVER_STARTED, sid + ": Getting a task.");

			Task task = Task.fromCommand(response);
			currentCalculationServer = ServerPool.getInstance().getFreeServer(task.taskType);
			redirectLog4jTo(currentCalculationServer.out);

			currentCalculationServer.setStatus(serverInfo.status = "Starting the task");
			calculationThread = new CalculationThread(task);
			calculationThread.start();

			// Notify ServerRunner that we are starting a task
			notifyRunner(ExchangeCommands.MULTISERVER_STARTED, sid + ": Started task " + task + " (" + task.getDataSize() % 1024 + "KB)");
			clean = false;
		}

		return false;
	}

	/**
	 * Append all Log4J events into the output stream
	 * 
	 * @param server
	 */
	public static void redirectLog4jTo(final PrintWriter targetWriter) {
		if (targetWriter == null)
			throw new RuntimeException("Cant redirect Log4J to a null PrintWriter");
		// FIMXE: Consider only removing the console appender and our own appenders. Allow other appenders to stay

		Logger logger = LogManager.getRootLogger();
		Configurator.setLevel(logger.getName(), Level.ERROR);

		/*

		LogManager.getRootLogger().removeAllAppenders();
		LogManager.getRootLogger().addAppender(new AppenderSkeleton() {

			@Override
			public boolean requiresLayout() {
				return false;
			}

			@Override
			public void close() {
			}

			@Override
			protected void append(LoggingEvent e) {
				targetWriter.println(e.getMessage());
			}
		});

		 */
	}

	public void taskCalculationMode() throws ClassNotFoundException, MalformedURLException, IOException, InterruptedException {

		// Notify the runner that MultiServer is still alive
		if (shoudRestart()) // corresponds to traditional mode -- server finishes after calculation
			notifyRunner(ExchangeCommands.MULTISERVER_PING, null);

		// No space left at all.. we have to stop the task
		if (!isEnoughDiskSpace(CRITICAL_DISK_SPACE)) {
			onTaskFinished(QSPRConstants.TASK_TEMPORAL_FAILURE + " Out of disk space. ");
			notifyRunner(ExchangeCommands.MULTISERVER_TERMINATE, "Out of disk space.. ");
		}

		serverInfo.cpuUsage = cpuMeasurementTool.getLoads()[2];

		// Task is running, report status to metaserver
		if (currentCalculationServer.status != null)
			serverInfo.status = currentCalculationServer.status;
		calculationThread.task.setDetailedStatus(serverInfo.status);

		Command response = transport.executeCommand(new Command(Command.CS_QUERYINFO, serverInfo).sid(sid));

		if (response == null)
			return;

		if(((Integer) response.data) == calculationThread.task.id)
			response.id = Command.MS_UNKNOWN_TASK;

		switch (response.id) {
		case Command.MS_UNKNOWN_TASK: 
			out.println("Metaserver does not know task: " + calculationThread.task + " -- it was deleted ");
			notifyRunner(ExchangeCommands.MULTISERVER_RESTART, "Killing task by request from server.");
			break;

		case Command.MS_STOP:
			currentCalculationServer.stopCalculations();
			break;

		case Command.MS_GET_LOGS:
			int debugLevel = calculationThread.task.debug;
			calculationThread.task.debug = DebugLevel.ALL; // Forcre logs extraction
			currentCalculationServer.takeLogsSnapshot(calculationThread.task.id);
			sendDebugFiles(calculationThread.task);
			calculationThread.task.debug = debugLevel;
			break;

		default:
			out.println("Ignoring command: " + response.id);
		}

	}

	// all communications are done through stderr,
	// to minimize problem with writing several messages simultaneously
	static private void notifyRunner(String command, String message) throws IOException {
		if (message != null && message.length() > 0) {
			System.err.println(message);
			out.println(message);
		}

		System.err.println("\n" + command);
		System.err.flush();

		if (!command.equals(ExchangeCommands.MULTISERVER_PING))
			out.println("\n" + command);

		// we should terminate or restart
		if (ExchangeCommands.MULTISERVER_RESTART.equals(command) || ExchangeCommands.MULTISERVER_TERMINATE.equals(command) ||
				ExchangeCommands.MULTISERVER_UPDATEREQUIRED.equals(command)
				)
			stopCalculationTask();
	}

	static private void stopCalculationTask(){

		if (calculationThread != null && calculationThread.isAlive()) {
			out.println("Interrupting task: " + calculationThread.task);
			calculationThread.interrupt();
		}

		ProcessUtils.killOffspring(null);

		if(Command.LOCAL.equals(version) || Command.FIXED.equals(version)) {
			out.println("Exiting with status 0");
			System.exit(0); // forcing to quit
		}else
			try {
				Thread.sleep(10000); // wait for 10 seconds and quit
				out.println("Not yet killed: harakiri with status 0");
				System.exit(0); // 
			} catch (InterruptedException e) {
			}
	}

	private void onTaskFinished(String message) throws MalformedURLException, IOException, ClassNotFoundException {
		// Task finished, post result back to metaserver
		// it may not work if metaserver is down (IOException)

		serverInfo.status = "Task finished";
		Task task = calculationThread.task;

		if (!calculationThread.finishedGracefully) 
		{
			task.setError((message != null ? message : "SEVERE: task failed presumably out of memory on ")
					+ ((serverInfo == null) ? "null" : serverInfo.name));
			out.println(task.getDetailedStatus());
		}

		sendDebugFiles(task);

		out.println("Task " + task + " finished. Status=" + task.status + "/" + task.getDetailedStatus() + "\nReporting to the metaserver...");
		task.peakMemoryUsage = MemoryUtils.getPeakMemoryUsage();
		task.clearData();
		transport.executeCommand(new Command(Command.CS_TASK_CALCULATED, task).sid(sid));
		serverInfo.status = "Finished task";
		calculationThread = null;
		out.println("Task " + task + " has been successfully reported to metaserver.");
	}

	private void sendDebugFiles(Task task)
	{
		if (task == null)
		{
			task = new Task();
			task.status = "error";
			task.debug = DebugLevel.ALL;
			task.taskType = "No task";
		}

		boolean sendLogsToMongodb = task.debug > DebugLevel.NONE && ("error".equals(task.status) || task.debug >= DebugLevel.ALL);
		File f = new File(currentCalculationServer.getDebugDirectory(task.id));
		if (f.exists())
		{
			File[] files = f.listFiles();
			for (File file : files) 
			{
				if (sendLogsToMongodb)
				{
					try
					{
						byte[] data = new byte[(int)file.length()];
						FileInputStream fis = new FileInputStream(file);
						fis.read(data);
						fis.close();
						Map<String, Object> map = new HashMap<String, Object>();
						map.put("filename", file.getName());
						map.put("task_id", task.id);
						map.put("parent_task_id", task.parentTaskId);
						map.put("server_address", serverInfo.ipAddress);
						map.put("server_name", serverInfo.name);
						map.put("server_configuration", serverInfo.configurationXml);
						map.put("task_type", task.taskType);
						NoSqlTransport.putDataSafely(data, "debug", "files", map);
						file.delete();
					}
					catch (Exception e)
					{
						e.printStackTrace(out);
					}
				}
			}
			f.delete();
		}
	}

	private boolean isEnoughDiskSpace(long limit) {
		long space = new File(workingDirectory).getUsableSpace();
		serverInfo.diskSpaceLeft = space / 1024 / 1024; // in MB
		boolean enough = space >= limit;

		if (!enough)
			out.println("SEVERE: There is not enough disk space to continue work! We have only " + space + " bytes");

		return enough;
	}

	private void configureSupportedServers(Applications applications, ServerRunnerConfiguration configuration) throws InstantiationException,
	IllegalAccessException, ClassNotFoundException, IOException {
		new CalculationServerConfigurator(serverInfo, out).configureServers(applications, configuration, workingDirectory);

		if (currentCalculationServer == null)
			currentCalculationServer = ServerPool.getInstance().servers.get(0);

		runTestsIfNecessary();
		serverInfo.setSupportedTasks(configuration.getSupportedTaskTypes(), configuration.currentVersion);
	}

	private void runTestsIfNecessary() throws IOException {
		if (skipTests)
		{
			out.println("Skipping tests as requested");
			return;
		}
		File fileWithFailureList = new File(workingDirectory + ExchangeCommands.MULTISERVER_TESTS_FILE_NAME);
		if (fileWithFailureList.exists())
			lastTestsRun = fileWithFailureList.lastModified();
		else
			out.println("No saved info about test failures. The tests will run now.");
		if (fileWithFailureList.exists() && (Calendar.getInstance().getTimeInMillis() - lastTestsRun < TEST_FAILED_RERUN_TIME)) // we use, of course, the least
			// time to rerun the test
		{
			// Tests have been already executed, and they are fresh (less than 1 hour old)
			// No need to run them again now

			BufferedReader reader = new BufferedReader(new FileReader(fileWithFailureList));
			String st;
			while ((st = reader.readLine()) != null)
				if (!st.equals("")) {
					serverInfo.addFailure(st);
					ServerPool.getInstance().disableTaskType(st);
				}
			reader.close();
			if (serverInfo.failures.size() > 0)
				out.println("Tests, failed on previous run: " + serverInfo.failures);
			else
				serverInfo.failures = null;
		} else {
			runTests();
			lastTestsRun = Calendar.getInstance().getTimeInMillis();

			// Write tests results to a file
			if (fileWithFailureList.exists())
				fileWithFailureList.delete();
			fileWithFailureList.createNewFile();
			BufferedWriter writer = new BufferedWriter(new FileWriter(new File(workingDirectory + ExchangeCommands.MULTISERVER_TESTS_FILE_NAME)));
			if (serverInfo.failures != null)
				for (String failure : serverInfo.failures) {
					writer.write(failure + "\n");
				}
			writer.flush();
			writer.close();
		}

		testsExpieryTime = serverInfo.failures == null ? TEST_PASSED_RERUN_TIME : TEST_FAILED_RERUN_TIME;
	}

	public void runTests() {
		out.println("Running tests...");
		JUnitCore jUnit = new JUnitCore();
		jUnit.addListener(new MyRunListener());
		if (generateTestReport) {
			try {
				jUnit.addListener(new ReportGeneratorTestListener(new File(workingDirectory + "/TEST-junit-report.xml")));
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		Request classes = null;
		try {
			classes = Request.classes(MainTest.class, Class.forName("qspr.metaserver.tests.EStateTest"), DescriptorServerTest.class, Class.forName("qspr.metaserver.tests.ScaffoldGeneratorTest"));
		} catch (Throwable exp) {
			classes = Request.classes(MainTest.class, DescriptorServerTest.class);
		}
		Request request = classes.filterWith(new Filter() {
			@Override
			public String describe() {
				return "Supported server tests";
			}

			@SuppressWarnings({ "unchecked", "rawtypes" })
			@Override
			public boolean shouldRun(Description description) {
				try {
					Class clazz;
					String name = description.getDisplayName(); 
					if(!name.contains("qspr.metaserver.tests.MainTest")) {
						out.println("skip: " + name);
						return false;
					}
					if (description.getDisplayName().contains("(")) {
						StringTokenizer tokenizer = new StringTokenizer(description.getDisplayName(), "(");
						String methodName = tokenizer.nextToken();
						String className = tokenizer.nextToken(")").substring(1);
						clazz = Class.forName(className);
						Method m = clazz.getMethod(methodName);
						return shouldRun(m);
					} else {
						clazz = Class.forName(description.getDisplayName());
						return shouldRun(clazz);
					}
				} catch (Exception e) {
					e.printStackTrace();
					return false;
				}
			}

			@SuppressWarnings("rawtypes")
			private boolean shouldRun(Class c) {
				Method[] methods = c.getDeclaredMethods();
				for (Method method : methods) {
					if (shouldRun(method))
						return true;
				}
				return false;
			}

			private boolean shouldRun(Method m) {
				TaskTest taskTest = (TaskTest) m.getAnnotation(TaskTest.class);
				if (taskTest != null)
					return ServerPool.getInstance().getFreeServer(taskTest.value()) != null;
				return false;
			}
		});
		jUnit.run(request);
	}

	class CalculationThread extends Thread {
		Task task;
		boolean finishedGracefully = false;
		long startedTime;

		public void run() {
			startedTime = Calendar.getInstance().getTimeInMillis();
			redirectLog4jTo(currentCalculationServer.out);
			currentCalculationServer.setStatus("Task started");
			currentCalculationServer.calculateWrapper(task);

			// If there had been some really bad uncaught runtime exception, like Java Heap Space,
			// this flag would stay false so that we could recognize failure / Midnighter
			finishedGracefully = true;

			redirectLog4jTo(out);
			// Notify the sleeping main thread, that the task is finished
			synchronized (this) {
				notifyAll();
			}

		}

		public long getRunningTimeSeconds() {
			return (Calendar.getInstance().getTimeInMillis() - startedTime) / 1000;
		}

		public CalculationThread(Task _task) {
			task = _task;
		}
	}

	class MyRunListener extends RunListener {
		private Map<String, Long> startTimes = new HashMap<String, Long>();

		public void testFailure(Failure failure) {
			// System.out.println(failure.getTestHeader()+" FAILED");
			// failure.getException().printStackTrace();
		}

		@Override
		public void testStarted(Description description) {
			out.println("Running test " + description.getDisplayName());
			startTimes.put(description.getDisplayName(), Calendar.getInstance().getTimeInMillis());
		}

		@Override
		public void testFinished(Description description) {
			Long timeStarted = startTimes.get(description.getDisplayName());
			if (timeStarted != null)
				out.println("Finished test " + description.getDisplayName() + " in  " + (Calendar.getInstance().getTimeInMillis() - timeStarted) + "ms.");
			else
				out.println("Finished test " + description.getDisplayName());
		}

		public void testRunFinished(Result result) throws Exception {
			out.println("Tests results:");
			for (Failure failure : result.getFailures()) {
				out.println(failure.getDescription().getDisplayName() + " FAILED - " + failure.getMessage());
				// if (failure.getDescription().getDisplayName().contains("Error"))
				failure.getException().printStackTrace(out);
			}
			out.println(result.getFailureCount() + " failures out of " + result.getRunCount() + " tests");

			// serverInfo.failures = new HashSet<String>();
			for (Failure failure : result.getFailures()) {
				String taskType = ((TaskTest) failure.getDescription().getAnnotation(TaskTest.class)).value();
				out.println(taskType);
				serverInfo.addFailure(taskType);

				// Now disable the server(s), that did not pass test
				ServerPool.getInstance().disableTaskType(taskType);
			}

			if (serverInfo.failures.size() == 0)
				serverInfo.failures = null;
			else
				out.println("List of disabled servers: " + serverInfo.failures);
		}
	}

	public static String getAutomaticSID(String workingDirectory) {
		String hostName = "";
		try {
			InetAddress localHost = InetAddress.getLocalHost();
			hostName = localHost.getHostName();
		} catch (UnknownHostException e) {
			// e.printStackTrace();
			hostName = e.getMessage();
		}
		String[] hostNameParts = hostName.split("\\.");
		String[] workingDirectoryParts = workingDirectory.split("\\/");
		String lastDirectoryPart = workingDirectoryParts[workingDirectoryParts.length - 1].replace("server", "");

		return hostNameParts[0] + "." + lastDirectoryPart;
	}

	private String getFileAsString(String file) throws IOException {
		BufferedReader input = new BufferedReader(new FileReader(file));
		StringBuffer sb = new StringBuffer();
		String line;
		while ((line = input.readLine()) != null)
			sb.append(line + '\n');
		input.close();
		return sb.toString();
	}
}

// Dublicates messages to System.out / Midnighter
class TunnelPrintWriter extends PrintWriter {
	String prefix;
	boolean timestamp = false;
	boolean silent = false;

	public PrintStream alternativeWriter = System.out;

	public TunnelPrintWriter(OutputStream out, String prefix) {
		super(out, true);
		if (prefix != null)
			this.prefix = "[" + prefix + "]";
		else
			this.prefix = "";
	}

	public TunnelPrintWriter(OutputStream out, String prefix, boolean timestamp, boolean silent) {
		super(out, true);
		if (prefix != null)
			this.prefix = "[" + prefix + "]";
		else
			this.prefix = "";
		this.timestamp = timestamp;
		this.silent = silent;
	}

	public TunnelPrintWriter setAlternativeWriter(PrintStream out)
	{
		this.alternativeWriter = out;
		return this;
	}

	public void println() {
		super.println();
		if (!silent)
			alternativeWriter.println();
	}

	public void print(String line) {

		String stamp = "";
		if (timestamp) {
			Calendar cal = Calendar.getInstance();
			SimpleDateFormat sdf = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.SSS");
			stamp = "[" + sdf.format(cal.getTime()) + "]";
		}
		super.print(stamp + " " + line);
		if (!silent)
			alternativeWriter.print(prefix + stamp + " " + line);
	}
}
