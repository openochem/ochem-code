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


import java.io.File;
import java.io.FilenameFilter;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.List;

import qspr.metaserver.protocol.Task;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.FileUtils;
import com.eadmet.utils.MemoryUtils;
import com.eadmet.utils.OCHEMUtils;

/**
 * Defines an abstract calculation servers that accepts a task (Task.java) and updates it with a calculated result TODO: Many things should be removed from this
 * class and moved down the hierarchy or to utility classes
 * 
 * @author midnighter
 * 
 */
public abstract class CalculationServer {

	public final String TOOLS ="tools"; // should be the same as in QSPRConstants
	static public final Integer NO_GPU = -1; 
	public final static String STOP = "stop";
	/**
	 * Indicates default task type for given server. Can be dynamically reconfigured using TASK in appplication.xml file
	 */
	public String supportedTaskType = null;
	/**
	 * Current status of server -- message that we see in metaserver
	 */
	public String status = "";
	public boolean busy = false;
	public static final String ERROR = "error";
	public PrintWriter out = new PrintWriter(System.out, true);
	public PrintWriter err = new PrintWriter(System.err, true);
	public static final String DATE_FORMAT_NOW = "dd-MM-yyyy HH:mm:ss";
	public static final String TEST_DIRECTORY = "test";
	/**
	 * The root directory of calculation server
	 */
	public String workingDirectory;
	/**
	 * The subdirectory with external binaries (those that are provided separately and can be configured by the user) 
	 */
	public String javaHome;
	public String javaClassPath;
	public Integer gpuCard = NO_GPU;
	public Integer availableMemory;

	Long start = null;

	private Double calculationProgress = null;

	public static boolean ENABLE_LICENSES = false; // Do not use any licenses
	protected String licensedModule = null;

	protected Event<String> statusChange = new Event<String>(this);

	/**
	 * This is directory in which calculations will run For server which do not use directory, it is a place in which one can store debugging information
	 */
	private String runsDirectory;

	/**
	 * The directory to pile-up all debug information, such as marshalled datatables for workflows, archives of particular "runs" and "output" directories, etc
	 */
	private String debugsDirectory;

	public CalculationServer() {
	}

	abstract protected void calculate(Task task) throws Exception;

	//	abstract protected boolean isInternal();

	protected int getMemoryForExecutable(){
		int memory = MemoryUtils.getCurrentMemoryFree()*3/5;
		return memory > availableMemory ? availableMemory : availableMemory/2; // at least half of the memory we will try to allocate
	}


	public String getToolsDirectory(){

		return workingDirectory + File.separator +TOOLS;
	}

	public void setStartTime() 
	{
		if (start == null) 
			start = System.currentTimeMillis();
	}

	public void resetStartTime() 
	{
		start = System.currentTimeMillis();
	}

	/**
	 * Set upper level directory that will be used to run programs (e.g., global place holders for all tasks) Directories for individual tasks will be created
	 * inside of it
	 * 
	 * @param d
	 */

	public void setRunsDirectory(String d) 
	{
		runsDirectory = d;
	}

	public void setDebugDirectory(String d)
	{
		debugsDirectory = d;
	}

	/**
	 * Provides upper level directory that is used to run program
	 * @return
	 */

	public String getRunsDirectory() 
	{
		return runsDirectory;
	}


	public String getDebugDirectory(Integer task_id)
	{
		String suffix;
		if (task_id == null)
			suffix = "global";
		else
			suffix = "task_"+task_id;
		File f = new File(debugsDirectory+suffix);
		if (!f.exists())
			f.mkdirs();
		return f.getAbsolutePath();
	}


	public void setParam(String name, String value) {
		out.println("Configuration option " + name + " = " + value);
		name = name.toUpperCase();
		if ("TASK".equals(name))
			supportedTaskType = value;
	}


	public void setStatus(String status) 
	{
		if (status == null) return;
		String fullStatus = status;
		if(status.length() > Task.MAX_STATUS_LENGTH)status = status.substring(0, Task.MAX_STATUS_LENGTH);

		if (!this.status.equals(status)) 
		{
			out.println(" [status] " + fullStatus);
			this.status = status;
			statusChange.fire();
		}
	}

	protected String getStatus() {
		return status;
	}

	/**
	 * Calculates a task and takes care of for errors and statistics
	 * @param task
	 */
	public final void calculateWrapper(Task task) {
		busy = true;
		long time = Calendar.getInstance().getTimeInMillis();
		out.println();
		out.println("Starting task " + task.id);
		out.println(MemoryUtils.memorySummary());
		try {

			calculate(task);
			task.status = "ready";
			task.setDetailedStatus("Finished");
			System.gc();
			out.println(MemoryUtils.memorySummary());
		} catch (InterruptedException e) {
			task.status = Task.KILLED;
			setStatus(task.setDetailedStatus("Task " + task + " has been killed"));
			out.println("Killing task " + task + ", because running thread has been interrupted");
		} catch (Exception e) {
			task.status = "error";
			if (e instanceof UserFriendlyException)
				task.setDetailedStatus(e.getMessage());
			else
				task.setDetailedStatus(OCHEMUtils.exceptionToString(e));
			setStatus(task.getDetailedStatus());
			out.println("Error in task " + task + ":\n" + task.getDetailedStatus());
			out.println("The exception that caused the failure: ");
			e.printStackTrace(out);
		} catch (Throwable e) 
		{
			// Something really bad happened here
			out.println("SEVERE: Something really bad has happened, an exception that could not be normally caught. ");
			setStatus(task.setDetailedStatus(e.getMessage()));
			System.gc();
			task.status = "error";
			task.setDetailedStatus(e.getMessage());
			e.printStackTrace(out);
			throw new AssertionError(e);
		} finally 
		{
			long interval = (Calendar.getInstance().getTimeInMillis() - time) / 1000;
			out.println("Finished task " + task.id + " in " + interval + " sec");
			busy = false;

			if (task.shouldKeepLogs())
			{
				takeLogsSnapshot(task.id);
			}

		}
	}

	public void takeLogsSnapshot(Integer task_id)
	{
		try
		{
			FileUtils.zip(new File(workingDirectory + "/output/"), new File(getDebugDirectory(task_id)+"/logs.zip"));
			if (task_id != null && new File(getTaskRunDirectory(task_id)).exists())
				FileUtils.zip(new File(getTaskRunDirectory(task_id)), new File(getDebugDirectory(task_id)+"/runs.zip"));

			// FIXME: Check deletion
		} catch (Exception e)
		{
			e.printStackTrace(out);
		}
	}

	/**
	 * 
	 */
	public final String getTaskRunDirectory(Integer taskId)
	{
		return runsDirectory
				+ File.separator
				+ supportedTaskType
				+ "_"
				+ (taskId == null ? TEST_DIRECTORY : taskId);
	}


	public static String now() {
		Calendar cal = Calendar.getInstance();
		SimpleDateFormat sdf = new SimpleDateFormat(DATE_FORMAT_NOW);
		return sdf.format(cal.getTime());
	}

	protected static MemoryPoolMXBean getPool() {
		List<MemoryPoolMXBean> pools = ManagementFactory.getMemoryPoolMXBeans();
		MemoryPoolMXBean result = null;
		for (MemoryPoolMXBean pool : pools) {
			if (pool.getName().contains("Old Gen")) {
				result = pool;
				break;
			}
		}
		return result;
	}

	public static long usedMemory() {
		MemoryPoolMXBean pool = getPool();
		return (pool == null) ? 0 : pool.getUsage().getUsed();
	}

	public static long maxUsedMemory() {
		MemoryPoolMXBean pool = getPool();
		return (pool == null) ? 0 : pool.getPeakUsage().getUsed();
	}

	public static void resetMaxUsedMemory() {
		MemoryPoolMXBean pool = getPool();
		if (pool != null)
			pool.resetPeakUsage();
	}

	/*
	 * Sets percentage of a task completed
	 */
	public void setPercentageCompleted(double val) 
	{
		calculationProgress = (val == 0) ? 0.001 : val;
	}

	/**
	 * 
	 * @return time required to complete task as String
	 */
	public String getTimeToComplete() {

		if (start == null)
			setStartTime();
		double seconds = 1.0 * (System.currentTimeMillis() - start) / 1000.0;
		if (calculationProgress == null)
			return null;

		if (seconds > 0) {

			seconds *= (1. / calculationProgress - 1.); // total duration - time already spent

			double t = seconds / (60 * 60 * 24.);

			if (t > 1)
				return "" + (int) (t * 10) / 10. + " days";
			t = seconds / (60 * 60.);
			if (t > 1)
				return "" + (int) (t * 10) / 10. + " hours";
			t = seconds / (60.);
			if (t > 1)
				return "" + (int) (t * 10) / 10. + " minutes";
		}

		return " a couple of minutes";
	}

	public void stopCalculations() { // we add stop to all directories in a hope that it will be useful
		File file = new File(runsDirectory);
		String[] directories = file.list(new FilenameFilter() {
			@Override
			public boolean accept(File current, String name) {
				return new File(current, name).isDirectory();
			}
		});
		for(String directory: directories) {
			FileUtils.touch(new File(runsDirectory+ "/" + directory + "/" + STOP));
		}
	}

	public void scanStatus(String status) {}

}
