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

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.HttpURLConnection;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.Unmarshaller;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.configuration.JAXBContextFactory;
import qspr.metaserver.protocol.Command;
import qspr.metaserver.protocol.Release;
import qspr.metaserver.transport.CSTransport;
import qspr.metaserver.util.AliasedDirectory;
import qspr.util.OverloadControl;
import qspr.workflow.utils.ProcessUtils;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.OSType;
import com.eadmet.utils.exe.OutReader;
import com.eadmet.utils.exe.ProperProcessClosing;

public class ServerRunner
{
	private static transient final Logger logger = LogManager.getLogger(ServerRunner.class);
	static final long CRITICAL_AVAILABLE_DISK_SPACE = 200 * 1024 * 1024; // 200 MB
	static final long UPDATE_TIME_REQUEST = 100; // seconds, considering that update itself requires about 3-5 minutes, we do not need to request it each second
	/**
	 * Maximum allowed time since last ping command from the MultiServer
	 */
	public static final String DATE_FORMAT_NOW = "dd-MM-yyyy HH:mm:ss";
	public static final String VERSIONXML = "version.xml";

	CSTransport transport;
	private ServerRunnerConfiguration configuration;
	Release localRelease;
	String serverHomeDirectory;
	File localConfigFile;

	JAXBContext jaxbContext;
	Unmarshaller unmarshaller;
	Marshaller marshaller;
	Process process;
	PrintWriter processPrinter;
	boolean runIt;
	private Thread mainThread;
	volatile long lastUpdateTime = 0;

	OverloadControl ocUpdateServer = new OverloadControl("update server", 5000, 60000);

	Release newRelease;
	/**
	 * Indicates that process is calculating a task now If false, server is just waiting for a task (idle mode)
	 */
	static volatile private boolean taskBeingCalculated = false;
	static volatile private long lastSignalFromServer;

	private long startupTime;
	private TunnelPrintWriter out;

	static final String METASERVER_IS_DOWN ="Server returned HTTP response code:";
	static final String JAVA_HEAP_SPACE ="java.lang.OutOfMemoryError: Java heap space";
	static final String TOOLS ="tools"; // should match same directory in CalculationServer

	public static void main(String[] args) throws Exception
	{
		ServerRunner m = new ServerRunner();
		try
		{
			m.initServer(args);
			m.out.println("ServerRunner is starting\n");
			m.run();
			m.marshaller.marshal(m.configuration, m.localConfigFile);
			logger.info("See you soon again.");
			System.exit(0);
		} catch (Exception e)
		{
			e.printStackTrace();
		}
	}

	private void initServer(String[] args) throws Exception
	{
		serverHomeDirectory = ".";
		mainThread = Thread.currentThread();
		startupTime = Calendar.getInstance().getTimeInMillis();
		localConfigFile = new File(serverHomeDirectory + "/" + VERSIONXML);

		// Do not run tests unless there is an update. No need. The test results will anyway expire
		// forceTestsRun(); // starting first time -- we should run tests

		if (args.length > 0)
			serverHomeDirectory = args[0];
		else
			serverHomeDirectory = new File(".").getCanonicalPath();
		if (!new File(serverHomeDirectory + File.separator + TOOLS).exists())
			throw new IOException("Invalid home directory: " + serverHomeDirectory + File.separator + TOOLS);

		out = new TunnelPrintWriter(new FileOutputStream(serverHomeDirectory + "/output/runner.out", true), null, true, false);
		MultiServer.redirectLog4jTo(out);
		ProcessUtils.out = out;

		out.println("Server Runner");
		UserCommand.printHelp(out);
		localRelease = new Release();

		jaxbContext = JAXBContextFactory.get("qspr.metaserver.serv");
		unmarshaller = jaxbContext.createUnmarshaller();
		marshaller = jaxbContext.createMarshaller();
		marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);

		ServerRunnerConfiguration.allApplications = getApplicationsFile();

		configuration = CalculationServerConfigurator.readConfigurationFile(jaxbContext, serverHomeDirectory);
		localRelease.version = configuration.currentVersion;

		out.println("Current version: " + configuration.currentVersion);

		if (configuration.metaserverURL != null && !configuration.metaserverURL.equals("")) {
			configuration.metaserverURL = CSTransport.addSeparator(configuration.metaserverURL);
			CSTransport.defaultServerURL = configuration.metaserverURL;
		}

		configuration.sid = MultiServer.getAutomaticSID(serverHomeDirectory);

		transport = new CSTransport();
		transport.serverURL = CSTransport.defaultServerURL + "/update";
		// Why ? If server is DEAD, not available, allow it some time to appear
		transport.deepSleepTime = 10 * 1000; // We use a sparing deep sleep interval here

		// Read input to execute user commands
		new OutReader(System.in, mainThread)
		{
			public void onLineRead(String line) throws Exception
			{
				executeUserCommand(line);
			}
		}.start();

		marshaller.marshal(configuration, localConfigFile);

		if (configuration.sleepTime < 1)
			throw new IOException("Sleep time should be at least 1 second");
	}

	/*	
	static public void changePermission(String path) throws IOException, InterruptedException{

		File file = new File(path + "/restart.sh");
		if(file.setExecutable(true))
			logger.info("Set executable status for "+file);
		else
			logger.info("Failed to set executable status for "+file);
	}

	static public void changePermissionOLd(String path) throws IOException, InterruptedException{

		Process proc = null;
		try
		{
			Runtime rt = Runtime.getRuntime();
			String[] commands = {"/bin/chmod"," u+x",path+"/*.sh"};
			proc = rt.exec(commands); // first make it executable
			proc.waitFor();
		} catch(Exception e){

		}
		finally
		{
			ProperProcessClosing.closeProcess(proc);
		}

	}
	 */	
	public void run() throws ClassNotFoundException, JAXBException, InterruptedException
	{

		runIt = true;

		int waitToKill = 0;

		// Main loop: check for updates, check if the server feels fine etc.
		// every second
		while (runIt)
		{
			// server is up and running (or waiting for a task)
			if (process != null)
				synchronized (mainThread)
				{
					mainThread.wait(configuration.sleepTime * 1000);
				}

			if (!taskBeingCalculated)
			{
				// Minimum lifetime reached?
				if (configuration.minimumLifetime != null && getUptime() > configuration.minimumLifetime)
				{
					out.println("Minimum lifetime reached. Kindly requesting the calculation server to restart...");
					if (isRunning(process))
					{
						requestTerminate();
						if(waitToKill++ * configuration.sleepTime > 3600){ 
							//Something strange -- server does not calculate task but process is still running (hanging)
							// We kill it after 1 hour
							out.println("NOT NORMAL: Maximum lifetime reached. Quitting.");
							killServer();
							break;
						}

						continue;
					}
					else
					{
						out.println("The server was not even running... We can safely quit now.");
						break;
					}
				}

				try
				{
					// If there is an update or server has not started -- start it!
					if (getUptime() - lastUpdateTime > UPDATE_TIME_REQUEST || process == null)
					{
						lastUpdateTime = getUptime();
						checkUpdate();
						if (process == null)
						{
							if (newRelease != null)
								performUpdate();
							startServer();
						}
					}
				} catch (IOException e)
				{
					// Check update or starting server has failed. Sleep, wait and retry again.
					e.printStackTrace(out);
					out.println("SEVERE: O-o it looks like that the update server is down. Waiting ...");
					lastSignalFromServer = Calendar.getInstance().getTimeInMillis(); // definitely in this case also CalculatIonServer cannot contact it!
				}
				continue; // We see that for this option we have finished work.
			}
			else
				waitToKill = 0; // task is calculated

			if (!isRunning(process))
			{
				// MultiServer was somehow killed or crashed for an unknown reason. Restart it.
				out.println("SEVERE: MultiServer was killed or crashed! Restaring it. Last signal was: "
						+ (Calendar.getInstance().getTimeInMillis() - lastSignalFromServer) / 1000 + " sec ago for process="+process);
				killServer(); // in case if nevertheless it is still alive :)
			}

			if (Calendar.getInstance().getTimeInMillis() - lastSignalFromServer > ExchangeCommands.SERVER_NO_ACTIVITY_MAX_TIME)
			{
				// The calculation server did not respond too long. That means a
				// crash or a hand of the main loop. Restart the server
				out.println("SEVERE: Server did not respond for > " + ExchangeCommands.SERVER_NO_ACTIVITY_MAX_TIME / 1000 + " secs. Restarting it...");
				killServer();
			}
		}

		out.println("ServerRunner exited normally. Bye-bye.");
	}

	public void executeUserCommand(String command) throws Exception
	{
		command = command.toUpperCase();
		String parts[] = command.split(" ", 2);
		UserCommand uCommand = null;
		try
		{
			uCommand = UserCommand.valueOf(parts[0]);
		} catch (IllegalArgumentException e)
		{
			out.println("Unrecognized command " + parts[0]);
			return;
		}

		switch (uCommand)
		{
		case HELP:
			out.println("Current version: " + configuration.currentVersion);
			UserCommand.printHelp(out);
			return;
		case RESTART:
			killServer();
			startServer();
			return;
		case RETEST:
			killServer();
			forceTestsRun();
			startServer();
			return;
		case QUIT:
			killServer();
			runIt = false;
			System.exit(0);
			return;
		case UPDATE:
			out.println("Forcing update");
			lastUpdateTime = 0;
			localRelease.version = "";
			interruptWaiting();
			return;
		case PRIORITY:
			configuration.minimumPriority = Integer.valueOf(parts[1]);
			marshaller.marshal(configuration, localConfigFile);
			out.println("Minimum priority has been set to " + configuration.minimumPriority);
			killServer();
			return;
		case LIST:
			Applications applications = getApplicationsFile();
			for (Application application : applications.applications)
			{
				out.print(application.name);
				if (configuration.isApplicationEnabled(application.name))
					out.print(" (TURNED ON)");
				out.println();
			}
			return;
		case ENABLE:
			if (configuration.enableApplication(parts[1]))
			{
				out.println("Application " + parts[1] + " has been enabled");
				marshaller.marshal(configuration, localConfigFile);
				killServer();
				forceTestsRun();
				startServer();
			}
			return;
		case DISABLE:
			if (configuration.disableApplication(parts[1]))
			{
				out.println("Application " + parts[1] + " has been disabled");
				marshaller.marshal(configuration, localConfigFile);
				killServer();
				startServer();
			}
			return;
		case FORWARD:
			messageToCalculationServer(parts[1]);

		}
	}

	public synchronized void executeCommandFromMultiServer(String command)
	{
		lastSignalFromServer = Calendar.getInstance().getTimeInMillis();

		if (ExchangeCommands.MULTISERVER_PING.indexOf(command) != -1 || command.contains(METASERVER_IS_DOWN))
			return;

		out.println("COMMAND: " + command);

		if(command.contains(JAVA_HEAP_SPACE)){
			out.println("Restarting the calculation server due to " + JAVA_HEAP_SPACE + " error");
			killServer();
			return;
		}

		if (ExchangeCommands.MULTISERVER_RESTART.indexOf(command) != -1)
		{
			out.println("Restarting the calculation server, as requested");
			killServer();
			return;
		}

		if (ExchangeCommands.MULTISERVER_RETEST.indexOf(command) != -1)
		{
			out.println("Restarting and retesting the calculation server, as requested");
			forceTestsRun(); 
			killServer();
			return;
		}

		if (ExchangeCommands.MULTISERVER_TERMINATE.indexOf(command) != -1)
		{
			out.println("Terminating server runner.");
			killServer();
			runIt = false;
			return;
		}

		if (ExchangeCommands.MULTISERVER_UPDATEREQUIRED.indexOf(command) != -1)
		{
			out.println("Multiserver said us that an update is required. Scheduling for an update ASAP.");
			// MultiServer is waiting for an update. It is in a sleeping mode and there is no need to wait with its killing.
			// Importantly -- server will not be started again until  new update will be installed
			killServer(); 
			return;
		}

		if (ExchangeCommands.MULTISERVER_STARTED.indexOf(command) != -1)
		{
			out.println("Calculation server has taken a task");
			taskBeingCalculated = true;
			return;
		}

		out.println("SEVERE: ERROR: not supported command:" + command);
	}

	private void interruptWaiting()
	{
		synchronized (mainThread)
		{
			mainThread.notifyAll();
		}
	}

	/**
	 * Before killing of the server (after update or since time has gone) we would like to notify it that we are ready to kill it (we request it restart)
	 * also verify that it did not take a new task meanwhile
	 * 
	 */
	private void requestRestart()
	{
		out.println("Requesting the multiserver to restart");
		messageToCalculationServer(ExchangeCommands.SERVER_RESTART_COMMAND);
	}

	private void requestTerminate()
	{
		out.println("Requesting the multiserver to terminate");
		messageToCalculationServer(ExchangeCommands.SERVER_TERMINATE_COMMAND);
	}

	private void messageToCalculationServer(String message)
	{
		if (process == null)
			return;
		if (processPrinter == null)
			processPrinter = new PrintWriter(process.getOutputStream(), true);
		processPrinter.println(message);
	}

	public void checkUpdate() throws IOException, ClassNotFoundException, JAXBException
	{

		out.println("Checking for new update.. (current: " + localRelease.version + ") using: " + transport.serverURL);

		if(transport.serverURL.startsWith(CSTransport.UNDEFINED))return;

		// First of all, do we have enough disk space?
		long space = new File(serverHomeDirectory).getUsableSpace();
		if (space < CRITICAL_AVAILABLE_DISK_SPACE)
		{
			out.println("SEVERE: Not enough disk space to download a release. " + space / (1024 * 1024) + "MB is available while "
					+ CRITICAL_AVAILABLE_DISK_SPACE / (1024 * 1024) + "MB is required");
			executeCommandFromMultiServer(ExchangeCommands.MULTISERVER_TERMINATE);
			return;
		}

		// Second, do we have a new release
		Command response = transport.executeCommand(new Command(0, newRelease != null ? new Release(newRelease.version) : localRelease).sid(configuration.sid));

		if (response != null)
		{
			newRelease = (Release) response.data;
			out.println("New update available ("+newRelease.version+" of "+newRelease.fileSize+" bytes)! Will download it at the next convenient time. ");
			out.println(MemoryUtils.memorySummary());
			// There is an update; we request server to terminate the server and will wait until it will finish
			requestRestart();
		}

	}

	private void getHttpToFile(String serverURL, File file) throws IOException
	{

		HttpURLConnection conn = null;		
		BufferedOutputStream bos = null;
		BufferedInputStream bis = null;

		IOException ee = null;

		try
		{
			bos = new BufferedOutputStream(new FileOutputStream(file));

			conn = CSTransport.getConnection(serverURL);

			conn.addRequestProperty("User-Agent", "OCHEM server"); // just to have this field not null
			conn.setReadTimeout(CSTransport.TRANSPORT_TIMEOUT); // 4 minutes
			conn.setConnectTimeout(CSTransport.TRANSPORT_TIMEOUT); // connection timeout
			conn.setRequestProperty("Connection", "close");

			bis = new BufferedInputStream(conn.getInputStream());

			byte[] buffer = new byte[4096];
			int numRead = 0;

			while ((numRead = bis.read(buffer)) > 0)
				bos.write(buffer, 0, numRead);

		}
		catch(IOException e){
			ee=e;
		}
		finally
		{
			conn.disconnect();
			if(bis!=null)bis.close();
			else
				logger.info("Exception was "+ee.getMessage());
			bos.close();
		}
	}

	synchronized public void performUpdate() throws IOException, InterruptedException
	{
		// Third, we kill running server, old release and install a new one
		try
		{
			if (configuration.currentVersion.equals(newRelease.version))
				return;

			File newReleaseArchive = new File(serverHomeDirectory + "/" + newRelease.version);
			out.println("Downloading new release "+newRelease.version);
			long time = System.nanoTime();
			getHttpToFile(CSTransport.defaultServerURL+"/update?" + configuration.sid, newReleaseArchive);
			out.println("Downloaded new release "+newRelease.version+". Total "+newReleaseArchive.length()+" bytes in "+(System.nanoTime() - time)*1.0/(1000*1000*1000)+" s");

			if(newReleaseArchive.length()==0)
				throw new IOException("new release failed to be downloaded -- will try again");

			if (newReleaseArchive.length() != newRelease.fileSize)
				throw new IOException(String.format("Corrupt release file %s: expected %d bytes and downloaded %d bytes", newRelease.version, newRelease.fileSize, newReleaseArchive.length()));

			out.println("Cleaning directory before the update ...");
			cleanDirectory();

			out.println("Applying update: " + newRelease.version);
			out.println("Unzipping using java ...");
			unzipFile(newReleaseArchive.getAbsolutePath());
			out.println("Unzip finished - update ready");
			configuration.currentVersion = localRelease.version = newRelease.version;

			newRelease = null;
			newReleaseArchive.delete();

			marshaller.marshal(configuration, localConfigFile);

			ServerRunnerConfiguration.allApplications = getApplicationsFile();

			// After update all tests should be re-run.
			forceTestsRun();

		} catch (Exception e)
		{
			e.printStackTrace(out);
			ocUpdateServer.relax(e);
			System.err.println(e.getMessage());
			configuration.currentVersion = localRelease.version = "error"; // we will need to fetch again
			throw new IOException("Update failed!");
		} finally
		{
			newRelease = null;
			System.gc();
		}

	}

	public synchronized void startServer() throws IOException
	{
		if (!runIt)
		{
			out.println("Calculation server will not be started, since the main loop is to be terminated");
			return;
		}

		if (process != null)
		{
			out.println("SEVERE: Server is already running");
			return;
		}

		logger.info("Starting the calculation server...");

		taskBeingCalculated = false;

		if (configuration.javaHome == null)
			configuration.javaHome = System.getProperty("java.home");

		int memory = configuration.memoryLimit > QSPRConstants.JAVA_MIN_MEMORY? configuration.memoryLimit : QSPRConstants.JAVA_MIN_MEMORY;

		String commands[] = {OSType.getPath(configuration.javaHome, "bin", "java"), "-cp",  // using Java as specified in cfg file; can be different from the one running ServerRunner
				OSType.getPath(serverHomeDirectory, "lib", "*"), "-Xmx" + memory +"m",
				"-XX:MinHeapFreeRatio=20", // an attempt to release memory aggressively by making it available to the native calculation servers
				"-XX:MaxHeapFreeRatio=40", // an attempt to release memory aggressively by making it available to the native calculation servers
				"qspr.metaserver.serv.MultiServer",serverHomeDirectory};

		process = Runtime.getRuntime().exec(commands);

		//out.println(commands);

		new OutReader(process.getErrorStream(), mainThread)
		{
			public void onLineRead(String line) throws Exception
			{
				// it does not matter which message is coming -- server is still
				// alive, if it can send it

				if (line.startsWith(ExchangeCommands.MULTISERVER_COMMAND_PREFIX) 
						|| line.contains(METASERVER_IS_DOWN) || line.contains(JAVA_HEAP_SPACE)) // special handling of errors will be done for them
				{
					executeCommandFromMultiServer(line);
				}

				if (!line.contains(ExchangeCommands.MULTISERVER_PING) && !line.isEmpty())
				{
					super.onLineRead("[" + ServerRunner.now() + "] " + line);
				}


			}
		}.setPreffix("[serv stderr:" + configuration.sid + "] ").start();

		// Buffering is required since some programs produce a lot of messages
		BufferedInputStream b = new BufferedInputStream(process.getInputStream());
		new OutReader(b, mainThread)
		{
			public void onLineRead(String line) throws Exception
			{
				// it does not matter which message is coming -- server is still
				// alive, if it can send it
				super.onLineRead(line);
			}
		}.setPreffix("[serv:" + configuration.sid + "] ").start();

		out.println("Server is running");

		if (processPrinter != null)
			processPrinter.close();
		processPrinter = null;

	}

	public synchronized void killServer()
	{
		taskBeingCalculated = false;
		if (process != null)
		{
			ProcessUtils.killOffspring(process);
			ProperProcessClosing.closeProcess(process);
			out.println("Server and its sub-processes have been killed");
		}
		else
			out.println("Killing ... but process was already null!");

		process = null;

		// after killing servers we may not have anymore time to work (including
		// time required to restart the server): thus finishing
		if (configuration.minimumLifetime != null && getUptime() > configuration.minimumLifetime * 8 / 10)
		{
			out.println("More than 80% of minimum lifetime has reached. Finishing calculations immediately ...");
			runIt = false;
		}

		interruptWaiting();

	}

	private boolean isRunning(Process process)
	{
		if (process == null)
			return false;

		boolean res = false;
		try
		{
			int exitvalue=process.exitValue();
			out.println("server exited with value="+exitvalue);
		} catch (IllegalThreadStateException e)
		{
			res = true;
		}
		return res;
	}

	private void unzipFile(String fileName) throws IOException
	{
		File file = new File(fileName);
		ZipFile zipFile = new ZipFile(fileName);
		Enumeration<? extends ZipEntry> entries = zipFile.entries();

		while (entries.hasMoreElements())
		{
			ZipEntry entry = (ZipEntry) entries.nextElement();

			if (entry.isDirectory())
			{
				// Assume directories are stored parents first then children.
				System.err.println("Extracting directory: " + entry.getName());
				// This is not robust, just for demonstration purposes.
				(new File(file.getParent() + "/" + entry.getName())).mkdir();
				continue;
			}

			System.err.println("Extracting file: " + entry.getName());
			copyInputStream(zipFile.getInputStream(entry), new BufferedOutputStream(new FileOutputStream(file.getParent() + "/" + entry.getName())));
		}

		zipFile.close();
	}

	public static String now()
	{
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat(DATE_FORMAT_NOW);
		return sdf.format(cal.getTime());
	}

	public void cleanDirectory() throws IOException, InterruptedException
	{
		// Delete all previous files and cleaning directory
		File[] files = new File(serverHomeDirectory).listFiles();
		for (File file : files)
			try
		{
				if (file.isDirectory())
				{
					out.println("Deleting directory: " + file.getName());
					AliasedDirectory.recursiveDelete(file);
				}
				else if (!file.getName().contains(VERSIONXML) && !file.getName().endsWith(".zip"))
				{
					out.println("Deleting " + file.getName());
					file.delete();
				}

		} catch (Exception e)
		{
			System.err.println("Exception to delete " + file.getName());
		}

	}

	public void forceTestsRun()
	{
		new File(serverHomeDirectory + ExchangeCommands.MULTISERVER_TESTS_FILE_NAME).delete();
	}

	private static final void copyInputStream(InputStream in, OutputStream out) throws IOException
	{
		byte[] buffer = new byte[1024];
		int len;

		while ((len = in.read(buffer)) >= 0)
			out.write(buffer, 0, len);

		in.close();
		out.close();
	}

	public long getUptime()
	{
		return (Calendar.getInstance().getTimeInMillis() - startupTime) / 1000;
	}

	private Applications getApplicationsFile() throws JAXBException
	{
		return (Applications) jaxbContext.createUnmarshaller().unmarshal(this.getClass().getClassLoader().getResourceAsStream("applications.xml"));
	}

}

class ReleaseFileFilter implements FilenameFilter
{
	public boolean accept(File dir, String name)
	{
		return (name.endsWith(".zip.tmp") || name.endsWith(".zip"));
	}
}

class Worker extends Thread
{
	private final Process process;
	private Integer exit;

	public Worker(Process process)
	{
		this.process = process;
	}

	public void run()
	{
		try
		{
			exit = process.waitFor();
		} catch (InterruptedException ignore)
		{
			return;
		}
	}

	public void waitFor(long timeout) throws InterruptedException, TimeoutException
	{
		start();
		join(timeout);
		if (exit == null)
		{
			process.destroy();
			throw new TimeoutException();
		}

	}
}

/**
 * Enumeration of supported command line commands
 * 
 * @author midnighter
 * 
 */
enum UserCommand
{
	QUIT("Quit the server completely"), RESTART("Restart the server without re-running tests"), RETEST("Restart the server with re-running tests"), UPDATE(
			"Force applying the new release"), LIST("List all available applications"), ENABLE(
					"[taskname] Enable task [taskname] and restart the server if necessary"), DISABLE(
							"[taskname] Disable task [taskname] and restart the server if necessary"), HELP("Display help"), PRIORITY(
									"Set the minimum priority for accepted tasks"),
	FORWARD("Forward a command to the downstream calculation server");

	public String description;

	private UserCommand(String description)
	{
		this.description = description;
	}

	public static void printHelp(PrintWriter out)
	{
		out.println();
		out.println("List of supported commands for the Server Runner");
		out.println("------------------------------------------------");
		for (UserCommand command : UserCommand.values())
			out.println(String.format("%-15s %s", command.name(), command.description));
	}
}

class TimeoutException extends Exception
{
	private static final long serialVersionUID = 1L;
};
