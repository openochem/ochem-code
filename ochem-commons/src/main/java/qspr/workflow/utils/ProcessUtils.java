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

package qspr.workflow.utils;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Map;
import java.util.Random;
import java.util.Timer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.eadmet.utils.FileUtils;
import com.eadmet.utils.OSType;
import com.eadmet.utils.exe.ProperProcessClosing;
import com.jezhumble.javasysmon.JavaSysMon;
import com.jezhumble.javasysmon.ProcessInfo;

public class ProcessUtils
{
	static JavaSysMon monitor = new JavaSysMon();

	private static transient final Logger logger = LogManager.getLogger(ProcessUtils.class);

	private static final String NSLOTS = "NSLOTS"; // number of slots available for server

	public static PrintWriter out = new PrintWriter(System.out);

	public final static String STDOUT ="stdout";

	public static void killByName(String aliasedPath) {
		logger.info("Killing all sub-processes by name " + aliasedPath);
		aliasedPath = aliasedPath.replaceAll("/", "").toLowerCase();
		logger.info("Killing all sub-processes by name " + aliasedPath);

		try {
			String line;
			Runtime rt = Runtime.getRuntime();
			Process p = rt.exec("ps -wax");
			BufferedReader input =
					new BufferedReader(new InputStreamReader(p.getInputStream()));
			while ((line = input.readLine()) != null) {
				String original = line;
				line = line.replaceAll("/", "").toLowerCase();
				if(line.contains(aliasedPath)) {
					String process[] = line.split("\\s+");
					int i = process[0].length()>0?0:1;
					logger.info("Killing  " + process[i] + " " + original);
					rt.exec("kill -9 "  + process[i]);
				}
			}
			input.close();
		} catch (Exception err) {
			err.printStackTrace();
		}		
	}

	// 
	/**
	 * It will kill children processes of the current process (not the process itself)
	 * Actually for Java 1.8 it will play mainly a role of waitFor
	 */
	public static void killOffspring(Process proc) 
	{
		logger.info("Killing all sub-processes of the current PID="+getCurrentPid());

		try{
			if(proc != null) proc.destroyForcibly(); // Java 8 implementation
		}catch(RuntimeException e){
			logger.info("proc.destroyForcibly() is not working ...");
		}catch(NoSuchMethodError e){
			logger.info("proc.destroyForcibly() is absent. Do we run Java 7?? ...");
		}	
		try{
			monitor.infanticide();
		}catch(RuntimeException e){
			logger.info("monitor.infanticide() is not working ...");
		}

		int counter = 0;
		try
		{
			while (monitor.processTree().find(getCurrentPid()).children().size() > 0)
			{
				counter++;
				logger.info("Waiting for child processes to finish, attempt "+counter);
				Thread.sleep(2000);
				if (counter > 10)
				{
					logger.error("Could not kill children in "+counter+" attempts." + proc == null? " and process was null. Really strange! ":"");
					break;
				}
			}
		} catch (Exception e)
		{
			logger.error("Thread sleep interrupted", e);
		}

		logger.info("Killing all sub-processes of the current PID="+getCurrentPid()+ " finished");

	}

	public static int getCurrentPid() 
	{
		return monitor.currentPid();
	}

	public static int getParentPid() 
	{
		ProcessInfo pr = getProcess(getCurrentPid());
		return pr == null ? null : pr.getParentPid();
	}

	public static boolean isProcessRunning(int pid)
	{
		return isProcessRunning(pid, null);
	}

	public static boolean isProcessRunning(int pid, Integer ppid)
	{
		ProcessInfo processes[] = monitor.processTable(); 
		for (ProcessInfo pr : processes)
		{
			if (pr.getPid() == pid && (ppid == null || ppid.equals(pr.getParentPid())))
			{
				logger.info("Found process ID " + pid + ", command " + pr.getCommand() + ", name " + pr.getName());
				return true;
			}
		}
		return false;
	}

	public static ProcessInfo getProcess(int pid) {
		ProcessInfo processes[] = monitor.processTable(); 
		for (ProcessInfo pr : processes)
		{
			if (pr.getPid() == pid)
				return pr;
		}
		return null;
	}

	public static String getComputerOwner()
	{
		return System.getProperty("user.name");
	}

	public static int getThreads(){
		Map<String, String> env = System.getenv();
		if(!env.containsKey(NSLOTS)) return 1;
		return Integer.parseInt(env.get(NSLOTS));
	}

	public static String runBinary(String[] commands) throws IOException, InterruptedException{
		return runBinary(commands, null, 1, null);
	}

	public static String runBinary(String[] commands,  String outputfile) throws IOException, InterruptedException{
		return runBinary(commands, outputfile, 1, null);
	}


	public static String runBinary(String[] commands,  String outputfile, int timeout_seconds) throws IOException, InterruptedException{
		return runBinary( commands, outputfile, timeout_seconds, null);
	} 

	public static void makeExecutable(String file) {
		try {
			File f = new File(file);
			if(f.exists() && !f.isDirectory()) {
				Runtime rt = Runtime.getRuntime();
				Process proc = rt.exec("chmod u+x "+file);
				proc.waitFor();
			}
		} catch (Exception e) {
		}
	}

	public static String runBinary(String[] commands, String outputfile, int timeout_seconds,  String[] environment) throws IOException 
	{
		Timer timer = null;
		Process proc = null;
		int exitVal = 0;
		boolean timeoutReached = false;

		// cleaning previous result file if it exists
		if(outputfile != null) {
			File f = new File(outputfile);
			if (f.exists())
				f.delete();
		}

		ExternalOutReader readOut = null, readErr = null;
		try
		{
			Runtime rt = Runtime.getRuntime();

			logger.info("executing: " + Arrays.toString(commands));

			if(!OSType.isWindows())
				makeExecutable(commands[0]);

			System.gc();

			proc = rt.exec(commands, environment);

			InterruptTimerTask interrupter = null;

			if (timeout_seconds > 0)
			{
				timer = new Timer(true);
				interrupter = new InterruptTimerTask(Thread.currentThread());
				timer.scheduleAtFixedRate(interrupter, timeout_seconds * 1000, timeout_seconds * 1000); // will execute regular check of interruption
			}

			readOut = (ExternalOutReader) new ExternalOutReader(new BufferedInputStream(proc.getInputStream()), Thread.currentThread(), null)
					.setOutputWriter(null)
					.setPreffix("");

			readOut.setExternalMonitor(interrupter);

			if (outputfile != null && STDOUT.equals(outputfile)) {
				Random r = new Random();
				outputfile = "/tmp/" + STDOUT+ "_"+r.nextInt();
				readOut.saveOutTo(outputfile);
			}

			readErr = (ExternalOutReader) new ExternalOutReader(new BufferedInputStream(proc.getErrorStream()), Thread.currentThread(), null)
					.setOutputWriter(null)
					.setPreffix("");

			readOut.start();
			readErr.start();

			readOut.join();
			readErr.join();

			// Are there any errors?
			exitVal = proc.waitFor();

			String mess = "";
			if(exitVal != 0 && mess.length() > 0)
				throw new IOException(mess);
		} catch (InterruptedException e)
		{
			String message = e == null || e.getMessage() == null?"":e.getMessage();
			logger.info("Timeout reached " + message);
			e.printStackTrace();
			timeoutReached = true;
			ProcessUtils.killOffspring(proc);
			proc.destroy();
		} catch (Exception e)
		{
			logger.info("error occurred while running the program: " + e.getMessage());
			throw new IOException(e.getMessage());
		} finally
		{

			if (readErr != null) 
				readErr.close();

			if (readOut != null) {
				readOut.flush();
				if (!timeoutReached) try{ // otherwise it hangs
					readOut.close();
				}catch(Exception e1) {}
			}

			if (timer != null)
			{
				timer.cancel();
				Thread.interrupted();
			}

			ProperProcessClosing.closeProcess(proc);
		}

		if(outputfile == null) return null;

		File outputFile = new File(outputfile);

		if (!outputFile.exists() || outputFile.length() == 0){
			String error = "External application of task crashed (non zero exit code " +exitVal+") and did not provide results";
			logger.info(error);
			throw new IOException(error);
		}

		return FileUtils.getFileAsString(outputFile.getPath());
	}

	public static void main(String[] args) throws Exception
	{
		String r[] = {"perl","/Users/itetko/a.pl"};
		String s = runBinary(r,STDOUT,10);

		System.out.println(s);

		/*		ProcessInfo processes[] = monitor.processTable(); 

		for (ProcessInfo pr:processes){
			System.out.println(pr);
			System.out.println(pr.getCommand());
			System.out.println(pr.getParentPid());

		}	 
		// A test
		logger.info(getComputerOwner());
		 */
	}

}
