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

package qspr.metaserver.util;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Timer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.dao.Various;
import qspr.metaserver.cs.WorkflowNodeServer;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.ExternalOutReader;
import qspr.workflow.utils.InterruptTimerTask;
import qspr.workflow.utils.ProcessUtils;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.CriticalException;
import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.FileUtils;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.OSType;
import com.eadmet.utils.exe.ProperProcessClosing;

/**
 * A helper class to run executable binaries.
 * 
 * @author itetko
 */
public class ExecutableRunner 
{
	protected static transient final Logger logger = LogManager.getLogger(ExecutableRunner.class);

	//	private String runsDirectory;

	WorkflowNodeServer parentServer;

	String BASHFILE = "bashCommand.sh";

	private String exefile = null;
	public final static String datain = "datain";
	public final static String dataout = "dataout";
	public final static String stdout = "stdout";

	public boolean logExeStdout = true;
	public boolean logExeSrderr = true;
	public String stdOutPrefix = "[exe stdout] ";
	public String stdErrPrefix = "[exe stderr] ";

	private AliasedDirectory aliasDirectory = null;

	public enum CONDA {RDKIT, MAP4, DEEPCHEM, CDDD, NULL};

	private String environmentConda(CONDA rdkit) throws IOException {

		switch(rdkit) {
		case DEEPCHEM: return "deepchem";
		case CDDD: return "cddd";
		case MAP4: return "map4";
		case RDKIT: return "my-rdkit-env";
		default:
			break;
		}

		throw new IOException("environment "+ rdkit +" not defined");
	}

	public ExecutableRunner(WorkflowNodeServer parentServer)
	{
		this.parentServer = parentServer;

		if(parentServer.params.containsKey(QSPRConstants.EXE))
			exefile = parentServer.params.get(QSPRConstants.EXE);
	}

	public void setLogStdOut(boolean val)
	{
		this.logExeStdout = val;
		parentServer.out.println("Logging exe stdout set to "+val);
	}

	public void setLogStdErr(boolean val)
	{
		this.logExeStdout = val;
		parentServer.out.println("Logging exe stderr set to "+val);
	}

	public AliasedDirectory getAliasDirectory()
	{
		return aliasDirectory;
	}


	/**
	 * Provides initialization of directories and possibly other common stuff
	 * required to use ServerAbstract This function should be called before
	 * using any other function.
	 * @throws InterruptedException 
	 */

	public void init() throws IOException, InterruptedException
	{
		Integer id = parentServer.currentTask.get().id;
		String taskRunsDirectory = parentServer.getTaskRunDirectory(id);

		if (taskRunsDirectory.startsWith("classpath:"))
		{
			parentServer.out.println("Skipping aliased directory initialization in classpath-home-girectory scenario");
			return;
		}
		aliasDirectory = new AliasedDirectory(exefile, parentServer.out, parentServer.getTaskRunDirectory(parentServer.currentTask.get().id));

		if(getExeFile() != null)
			parentServer.out.println("aliasedBinary: " + getExeFile());
		parentServer.out.println("aliasedDirectory: " + aliasDirectory.getAliasedPath());
	}

	/**
	 * Delete all temporal directories
	 * @throws InterruptedException 
	 */

	public void cleanup() throws IOException, InterruptedException
	{
		if (aliasDirectory == null)
			return;

		if (!parentServer.currentTask.get().shouldKeepLogs())
		{
			parentServer.out.println("delete aliasedDirectory: " + aliasDirectory.getAliasedPath());
			aliasDirectory.deleteAlias();
		}
	}

	public String getExeFile()
	{
		return exefile == null? null:getAliasedFileName(exefile);
	}

	public void setExeFile(String file) throws IOException
	{
		exefile = file;
		File f = new File(exefile);
		if (f.exists() && f.isFile())
			f.setExecutable(true);
		else
			throw new IOException("Provided executable file "+f.getAbsolutePath()+" does not exist or is not a file");
	}

	public String getAliasedFileName(String filename)
	{
		if (aliasDirectory == null)
			throw new RuntimeException("Directory to execute program is not specified. Did you call init function?");
		return aliasDirectory.getAliasedPath(filename);
	}

	public String getAliasedPath()
	{
		if (aliasDirectory == null)
			throw new RuntimeException("Directory to execute program is not specified. Did you call init function?");
		return aliasDirectory.getAliasedPath();
	}
	/**
	 * Creates file reader from the aliased filename
	 * 
	 * @param filename
	 * @return
	 * @throws Exception
	 */
	public BufferedReader getAliasedBufferedReader(String filename) throws IOException
	{
		FileReader fis = new FileReader(getAliasedFileName(filename));
		BufferedReader bis = new BufferedReader(fis,QSPRConstants.FILE_BUFFER);
		return bis;
	}

	/**
	 * Creates file writer from the aliased filename
	 * 
	 * @param filename
	 * @return
	 * @throws Exception
	 */
	public BufferedWriter getAliasedBufferedWriter(String filename) throws IOException
	{
		FileWriter fis = new FileWriter(getAliasedFileName(filename));
		BufferedWriter bis = new BufferedWriter(fis,QSPRConstants.FILE_BUFFER);
		return bis;
	}

	// Command to run the program
	// as well as file that will be used to check whether the program has
	// completed

	private String toString(String[] a)
	{
		String r = "";
		for (String aa : a)
		{
			r += aa + " ";
		}
		return r;
	}

	/**
	 * Executes commands provided by commands Uses outputfile to verify that
	 * calculations were executed successfully if outputfile name is equal to
	 * stdout, catches stdout and write it to this file
	 * 
	 * @param commands
	 * @param res
	 * @throws InterruptedException 
	 * @throws Exception
	 */
	public void executeBinary(String[] commands, String outputfile) throws IOException, InterruptedException
	{
		executeBinary(commands, outputfile, 0);
	}

	/**
	 * 
	 * @param commands
	 * @param outputfile - file that will be created as result of this commands
	 * @param timeout_seconds - timeout in seconds
	 * @throws IOException
	 * @throws InterruptedException
	 */

	public void executeBinary(String[] commands, String outputfile, int timeout_seconds) throws IOException, InterruptedException 
	{
		executeBinary(commands, null, outputfile, timeout_seconds, parentServer.out, parentServer.err);
	}

	/**
	 * @param commands
	 * @param outputfiles files that have to be eliminated before calculation; first one will be used to monitor that calculations are finished
	 * @param timeout_seconds is in seconds
	 * @param outputWriter
	 * @throws IOException
	 * @throws InterruptedException 
	 */
	public void executeBinary(String[] commands, String[] environment, String outputfiles, int timeout_seconds, PrintWriter outputWriter, PrintWriter errWriter) throws IOException, InterruptedException 
	{
		if (aliasDirectory == null || aliasDirectory.getAliasedPath() == null)
			throw new IOException("Directory to execute program is not specified. Did you call init function?");

		Timer timer = null;
		Process proc = null;
		int exitVal = 0;
		boolean timeoutReached = false;
		ProcessUtils.killByName(parentServer.getRunsDirectory());

		String files[] = outputfiles.split(","); 
		for(String file : files) {
			// cleaning previous result file if it exists
			File f = new File(getAliasedFileName(file));
			if (f.exists())
				f.delete();
		}
		outputfiles = files[0]; // only the first file will be checked that results were calcukateds

		ExternalOutReader readOut = null, readErr = null;
		try
		{
			parentServer.out.println("aliased dir: " + getAliasedPath());

			Runtime rt = Runtime.getRuntime();

			if(!OSType.isWindows())
				for(String command : commands) {
					if(command.contains(BASHFILE)){
						String s = FileUtils.getFileAsString(getAliasedFileName(command));
						for(String file:s.split("\\s+")) {
							ProcessUtils.makeExecutable(file);
						}	
					}
					else 
						ProcessUtils.makeExecutable(command);
				}

			parentServer.out.println("running commands: " + toString(commands));

			parentServer.setStartTime();

			System.gc();

			if(ProcessUtils.getThreads()>1) {
				environment = OCHEMUtils.appendString(environment,QSPRConstants.OMP_NUM_THREADS+"=" + ProcessUtils.getThreads());
				parentServer.out.println("Using parallel environment " +  environment[environment.length-1]);
			}

			proc = rt.exec(commands, environment, new File(getAliasedPath()));

			InterruptTimerTask interrupter = null;

			if (timeout_seconds > 0)
			{
				timer = new Timer(true);
				interrupter = new InterruptTimerTask(Thread.currentThread());
				timer.scheduleAtFixedRate(interrupter, timeout_seconds * 1000, timeout_seconds * 1000); // will execute regular check of interruption
			}

			readOut = (ExternalOutReader) new ExternalOutReader(new BufferedInputStream(proc.getInputStream()), Thread.currentThread(), parentServer)
					.setOutputWriter(logExeStdout ? outputWriter : null)
					.setPreffix(stdOutPrefix);

			readOut.setExternalMonitor(interrupter);

			if (stdout.equals(outputfiles))
				readOut.saveOutTo(getAliasedFileName(stdout));

			readErr = (ExternalOutReader) new ExternalOutReader(new BufferedInputStream(proc.getErrorStream()), Thread.currentThread(), parentServer)
					.setOutputWriter(logExeSrderr ? errWriter : null)
					.setPreffix(stdErrPrefix);

			readOut.start();
			readErr.start();

			readOut.join();
			readErr.join();

			// Are there any errors?
			exitVal = proc.waitFor();

			String mess = "";
			if (readOut.errorMessage != null || readErr.errorMessage != null){
				mess = readOut.errorMessage != null? readOut.errorMessage : "" + " " 
						+ readErr.errorMessage != null? readErr.errorMessage : "";
				readOut.errorMessage = readErr.errorMessage = null;
			}
			mess += readErr.getLastMessage();
			if(exitVal != 0 && mess.length() > 0)
				throw new IOException(mess);
		} catch (InterruptedException e)
		{
			String message = e == null || e.getMessage() == null?"":e.getMessage();
			logger.info("Timeout reached " + message);
			e.printStackTrace();
			parentServer.out.println(message);
			timeoutReached = true;
			parentServer.out.println("Timeout of " + timeout_seconds + " seconds was reached");
			ProcessUtils.killOffspring(proc);
			ProcessUtils.killByName(getAliasedPath());
			proc.destroy();
			parentServer.out.println("Killing external process -- finished");
		} catch (Exception e)
		{
			logger.info("No timeout!");
			logger.info("error occurred while running the program: " + e.getMessage());
			if(parentServer.isCritical(e.getMessage()))
				throw new CriticalException(e.getMessage());
			throw new UserFriendlyException(e.getMessage());
		} finally
		{

			if (readErr != null) {
				parentServer.out.println("closing readErr");
				readErr.close();
			}

			if (readOut != null) {
				readOut.flush();
				if(!timeoutReached) try{ // otherwise it hangs
					parentServer.out.println("closing readOut");
					readOut.close();
				}catch(Exception e1) {}
			}

			if (timer != null)
			{
				parentServer.out.println("stopping timer");
				timer.cancel();
				Thread.interrupted();
			}

			parentServer.out.println("closing process");
			ProperProcessClosing.closeProcess(proc);
		}

		File outputFile = new File(aliasDirectory.getAliasedPath(outputfiles));
		if (!outputFile.exists() || ( outputFile.length() == 0 && !outputFile.isDirectory()))
		{
			// only if output file is not provided we throw exception, since it
			// can be the reason why it is not available
			if (exitVal != 0)
				throw new UserFriendlyException("External application of task " +parentServer.supportedTaskType+" crashed (non zero exit code "
						+exitVal+") and did not provide results");
			throw new UserFriendlyException("External application of task " +parentServer.supportedTaskType+" did not provide results at " + outputFile);
		}

		if (exitVal != 0)
			logger.error(" ATTENTION: The program has presumably failed but still produced a file with results. Using it.");

		int counter = 0;
		while (!(outputFile.canRead()) && !(outputFile.canWrite()))
		{
			counter++;
			String message = outputFile.canRead() ? "output file is ready to read" : "output file is not readable";
			message += outputFile.canWrite() ? "output file is ready to write" : "output file is not writeable";

			parentServer.setStatus("Waiting for file to be ready: " + message + " " + counter + "s");
			if (counter > 1000)
				throw new IOException("reached maximun waiting time. program will terminate the model calculation. " + message
						+ " try later and still problem exist then please report this bug");
			try {
				Thread.sleep(100);
			} catch (InterruptedException e) {
			}
		}

	}

	/*
	 * Looks for the binary in pre-specified directories by name 
	 */

	public static String findExecutable(String executable)
	{
		String places[] = { "/opt/local/bin/", "/opt/conda/bin/", "/usr/local/bin/", "/usr/bin/"};

		for (String s : places)
			if (checkExecutable(s + executable))
				return s + executable;

		return executable; // not found, maybe somewhere else in path?
	}

	static boolean ignore = true;

	public String findWorkingPython(String[] commands, String results, int position, CONDA conda, int timeout) throws Exception{

		if(commands[position] != null) {
			executeBinaryBash(commands, results,conda,0);
			return commands[position];
		}

		String pythons [] = {};

		conda = conda == null ? CONDA.NULL : conda;

		if(OSType.isMac()) { // only for debugging...
			switch(conda) {
			case MAP4:
			case CDDD:
			case RDKIT:
				pythons = OCHEMUtils.append(pythons,OSType.getHome()+ QSPRConstants.TETKO_ANACONDA + "/envs/"+environmentConda(conda)+"/bin/pythonw"); 
			case DEEPCHEM:
				pythons = OCHEMUtils.append(pythons,OSType.getHome()+ QSPRConstants.TETKO_ANACONDA  + "/envs/"+environmentConda(conda)+"/bin/python"); break;
			default:
				//pythons = new String[] {"/opt/local/bin/python3.8","/opt/local/bin/python3.9"};
				pythons = new String[] {"/opt/local/bin/python3.10"};
			}
		}else {
			pythons = OCHEMUtils.append(pythons,"/opt/conda/bin/python");
			switch(conda) {
			case RDKIT:
				pythons = OCHEMUtils.append(pythons,"python");
				break;
			default:
				pythons = OCHEMUtils.append(pythons,"python");
				break;
			}

		}

		Exception ee = null;

		for(String s:pythons)try{
			String file;
			if((file = OCHEMUtils.findFile(s)) != null){
				File f = new File(file);
				commands[position] = OSType.isMac()?s:aliasDirectory.alias(new File(f.getParent()), new String[] {f.getName()});
				executeBinaryBash(commands, results, conda, timeout);
				return commands[position];
			}
		}catch(Exception e) {
			if(ee == null) ee = e;
		}
		if(ee == null)throw new IOException("required python could not be found");
		else
			throw ee;
	}

	static private boolean checkExecutable(String location) {
		File f = new File(location);
		return f.canExecute();
	}

	public void saveMolecules(DataTable dtMolecules, String filename,
			String formatToSave) throws IOException {
		saveMolecules(dtMolecules, filename, formatToSave, 0, dtMolecules.getRowsSize());
	}

	public void saveMolecules(DataTable dtMolecules, String fileName, String formatToSave, int rowFrom, int rowCount) throws IOException{
		BufferedWriter writer = getAliasedBufferedWriter(fileName);
		saveMolecules(dtMolecules, writer, formatToSave, rowFrom, rowCount);
		writer.close();
	}

	static public void saveMolecules(DataTable dtMolecules, BufferedWriter writer, String format, int rowFrom, int rowCount) throws IOException
	{

		String eL=OSType.endLine(),eL3=eL+eL+eL;

		if (format.equals(QSPRConstants.MOPIN))  //special handling of MOPAC format : it should be already in this format ...
		{
			if(rowCount-rowFrom>1)throw new IOException("Only 1 molecule is allowed to be saved for "+QSPRConstants.MOPIN+" format");
			writer.write((String) dtMolecules.getValue(QSPRConstants.MOPIN) + eL3);
		}
		else
			for (int i = rowFrom; i < rowFrom + rowCount && i < dtMolecules.getRowsSize(); i++)
			{
				// first, we always import the molecule
				// molecules can be in SDF or in SMILES

				String mol=(String) dtMolecules.getValue(i, 0);

				switch(format){

				case QSPRConstants.SDF: 
				case QSPRConstants.SDFH: 
				case QSPRConstants.SDFNOAROM_NOH: 
				case QSPRConstants.SDFAROM_BASIC_WITHH: 
				case QSPRConstants.SDFAROM_BASIC_NOH:
				case QSPRConstants.SDFAROM_GENERAL_WITHH:
				case QSPRConstants.SMILES_FORMAT:
				case QSPRConstants.SMILESH:
				case QSPRConstants.SMILESNOAROM:
				case QSPRConstants.SMILESNOSTEREO:
				case QSPRConstants.SMILESUniqueNoHAromatic:
					String smiles = Various.molecule.convertToFormat(mol,format);
					for(String line:smiles.split("\\r?\\n"))
						writer.write(line+eL);
					break;

				case QSPRConstants.SDFNOAROM_WITHH:
					String sdf = Various.molecule.addHydrogensAndRemoveRadicalsAndSMARTS(mol);
					for(String line:sdf.split("\\r?\\n"))
						writer.write(line+eL);
					break;

				case QSPRConstants.MOL2:
					String mol2 = Various.molecule.convertToFormat(mol,format);
					writer.write(mol2);
					break;

				case QSPRConstants.SMILESSRC:
					String converted = Various.molecule.convertToFormat(mol,QSPRConstants.SMILESNOSTEREO);
					String smile = converted.replaceAll("\\[N\\]", "(N)");
					if(!converted.equals(smile))System.out.println(converted+" "+smile);
					writer.write(smile + eL);
					break;

				case QSPRConstants.ASIS:
					for(String line:mol.split("\\r?\\n"))
						writer.write(line+eL);
					break;

				default: 
					throw new IOException("Format is not supported " + format);
				}

			}

		writer.close();
	}

	public String runPython(String [] commands, String file, CONDA conda, int timeout_seconds) throws Exception{

		int found = OCHEMUtils.fibdNullPosition(commands);
		if(found != -1) {
			parentServer.out.println("Starting calculations " + parentServer.supportedTaskType);
			findWorkingPython(commands, file, found, conda, timeout_seconds);
			parentServer.out.println("Found " + commands[found]);
			return commands[found];
		}
		else
			executeBinaryBash(commands, file, conda, timeout_seconds);

		return null;

	}

	public void executeBinaryBash(String[] commands, String outputFile, CONDA conda, int timeout_second) throws IOException, InterruptedException {
		executeBinaryBash(commands, null, outputFile, conda, timeout_second);
	}
	/*
	 * Run already prepared command including both 
	 */
	public void executeBinaryBash(String[] commands, String[] env, String outputFile, CONDA conda, int timeout_second) throws IOException, InterruptedException {

		BufferedWriter writer = getAliasedBufferedWriter(BASHFILE);
		if(conda != null && conda != CONDA.NULL) {
			if(OSType.isMac())
				writer.append("source " + OSType.getHome()+ QSPRConstants.TETKO_ANACONDA + "/bin/activate " + environmentConda(conda) +"; export KMP_DUPLICATE_LIB_OK=TRUE; ");
		}

		for(String s : commands) { 
			writer.append(s+" ");
			if(s.contains("python")) writer.append("-u ");
		}

		writer.close();
		String com[] = {"/bin/bash",BASHFILE};
		parentServer.out.println("running following commands with timeout " + timeout_second + " using bash file: " + toString(commands));
		executeBinary(com, env, outputFile, timeout_second, parentServer.out, parentServer.err);
	}

	public static void main(String args[]) {
		String smi="CC1=CC=CC=C1[N](O)=O";
		smi=smi.replaceAll("\\[N\\]", "(N)");
		System.out.println(smi);
	}
}