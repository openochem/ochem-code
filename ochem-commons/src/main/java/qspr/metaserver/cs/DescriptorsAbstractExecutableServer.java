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

package qspr.metaserver.cs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.LineNumberReader;
import java.io.PrintWriter;
import java.util.Calendar;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.util.ExecutableRunner;
import qspr.metaserver.util.ExecutableRunner.CONDA;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

public abstract class DescriptorsAbstractExecutableServer extends DescriptorsAbstractServer
{
	private DataTable dtResults;
	protected DataTable dtConditions; // required for calculation of some types of descriptors
	public final static String datain = "datain";
	public final static String dataout = "dataout";
	public final static String stdout = "stdout";
	protected String DELIMITER = "\\s";
	protected int startPosition = 1;

	protected ExecutableRunner exeRunner;

	void init() throws IOException, InterruptedException{
		exeRunner = new ExecutableRunner(this);
		exeRunner.init();
	}

	void cleanup() throws IOException, InterruptedException{
		exeRunner.cleanup();
		exeRunner=null;
	}

	protected String getAliasedPath() {
		return exeRunner.getAliasedPath();
	}

	protected String getAliasedFileName(String filename) {
		return exeRunner.getAliasedFileName(filename);
	}

	protected LineNumberReader getAliasedLineNumberReader(String filename) throws IOException {
		return new LineNumberReader(exeRunner.getAliasedBufferedReader(filename));
	}

	protected BufferedWriter getAliasedBufferedWriter(String filename) throws IOException {
		return exeRunner.getAliasedBufferedWriter(filename);
	}

	public BufferedReader getAliasedBufferedReader(String filename) throws IOException {
		return exeRunner.getAliasedBufferedReader(filename);
	}

	protected void runPython(String[] commands, String outputFile, CONDA rdkit, int timeout) throws Exception{
		exeRunner.runPython(commands, outputFile, rdkit, timeout);
	}

	protected void executeBinary(String[] commands, String outputFile) throws IOException, InterruptedException {
		exeRunner.executeBinary(commands, outputFile);
	}

	public void executeBinaryBash(String[] commands, String env[], String outputFile, int timeout_seconds) throws IOException, InterruptedException {
		exeRunner.executeBinaryBash(commands, env, outputFile, null, timeout_seconds);
	}

	public void executeBinaryBash(String[] commands, String outputFile, int timeout_seconds) throws IOException, InterruptedException {
		exeRunner.executeBinaryBash(commands, outputFile, null, timeout_seconds);
	}

	protected void executeBinary(String[] commands, String outputFile, int timeout_seconds) throws IOException, InterruptedException {
		exeRunner.executeBinary(commands, new String[] {}, outputFile, timeout_seconds, out, err);
	}

	protected void executeBinary(String[] commands, String env[], String outputFile, int timeout) throws IOException, InterruptedException {
		exeRunner.executeBinary(commands, env, outputFile, timeout, out, err);
	}

	protected void executeBinary(String[] commands, String outputFile, int timeout,
			PrintWriter out, PrintWriter err) throws IOException, InterruptedException {
		exeRunner.executeBinary(commands, new String[] {}, outputFile,timeout,out,err);
	}

	protected String getExeFile() {
		return  exeRunner.getExeFile();
	}

	protected DataTable getResults() throws IOException
	{
		if (dtResults == null)
			throw new IOException("dtResults is not yet initialised.");
		return dtResults;
	}

	@Override
	public final WorkflowNodeData calculateDescriptors(WorkflowNodeData task, DescriptorsAbstractConfiguration configuration) throws Exception
	{


		try {
			DataTable mols = task.ports.get(0);
			dtResults = new DataTable(isDataTableCompact());
			dtConditions = task.ports.size() > 1 ? task.ports.get(1) : null;
			DataTable dtDescriptors = null;

			init();
			dtDescriptors = calculateInBatch(mols, configuration, getBatchSize());
			fixNullMolecules(dtDescriptors, mols);
			if(!(isRunningTest() && dtDescriptors.getRowsNoErrorsSize() == 0)) // test failed, let us keep directory
				cleanup();

			dtDescriptors.id = dtDescriptors.id == null ? "descriptors" : dtDescriptors.id;
			return new WorkflowNodeData(dtDescriptors);
		}catch(Throwable e) {
			throw e;
		}finally {
			//Various.molecule = molecule;
		}
	}

	private void fixNullMolecules(DataTable dtDescriptors, DataTable dtMolecules) {
		for(int start = 0;  start < dtDescriptors.getRowsSize(); start++) 
			if( (dtMolecules.getValue(start, 0) == null ||  ((String)dtMolecules.getValue(start, 0)).length() == 0))
				dtDescriptors.getRow(start).setError(QSPRConstants.EMPTY_MOLECULE_ERROR);
	}
	/**
	 * Should be overwritten to provide information 
	 * about the batch size if processing of all molecules is long/crashes
	 * @return
	 */

	int getBatchSize()
	{
		return 0;
	}

	/**
	 * Should be overwritten to provide information 
	 * that data table is not compact
	 * @return
	 */
	protected boolean isDataTableCompact()
	{
		return true;
	}

	void saveMolecules(DataTable dtMolecules, String filename, String datatype,int from,int rowCount) throws IOException {
		exeRunner.saveMolecules(dtMolecules,filename,datatype,from,rowCount);		
	}

	void saveMolecules(DataTable dtMolecules, BufferedWriter writer, String datatype,int from,int rowCount) throws IOException {
		ExecutableRunner.saveMolecules(dtMolecules,writer,datatype,from,rowCount);		
	}

	/**
	 * Reads descriptors in standard format: 
	 * name desName1 desName2 desName3
	 * mol1 val11 val12 val13 
	 * mol2 val21 val22 val23 
	 * and so on 
	 * first column is ID of a molecule first row contain names of columns
	 * 
	 * @param filename
	 *            - file with descriptors
	 * @return DataTable with descriptors
	 * @throws Exception
	 *             if there is reading error
	 */

	protected DataTable readStandardOutputResults(DataTable dres,String filename, boolean forceAddDecriptorNames) throws IOException
	{
		BufferedReader outputReader = getAliasedBufferedReader(filename);

		String line;
		boolean first = true;
		while ((line = outputReader.readLine()) != null)
		{

			line = line.replaceAll("\"", "");

			String[] n = line.split(DELIMITER);

			//if(n.length > 108 && n[108].equals("GEIGEE_elcdloc")) // FIX
			//	n[108]="GEIGEE_nucdloc"; //FIX

			if (first)
			{
				first = false;

				if(forceAddDecriptorNames || dres.getColumnsSize() == 0)
					for (int i = startPosition; i < n.length; i++)
						dres.addColumn(n[i]);
				continue;
			}

			dres.addRow();

			if (line.toLowerCase().contains("error"))
				dres.getCurrentRow().setError(line.substring(line.indexOf("error") + "error".length() + 1));
			else
				if ((dres.getColumnsSize() ) != (n.length - startPosition))
					dres.getCurrentRow().setError("Expected: " + dres.getColumnsSize() + " get: " + n.length + " descriptor");
				else
					for (int i = startPosition; i < n.length; i++) {
						if(n[i].equalsIgnoreCase("NA"))
							n[i]="NaN";
						if(n[i].length()>0)dres.setValue(i - startPosition, Double.parseDouble(n[i]));
					}
		}

		outputReader.close();
		return dres;
	}

	/**
	 *  Calculates descriptors starting from Start to Stop
	 *  If fails, calculates them one by one
	 * @param dtMolecules
	 * @param dtConditions
	 * @param receivedConfiguration
	 * @param MolPerBatch
	 * @return
	 * @throws Exception
	 */

	private DataTable calculateInBatch(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration, int batch) throws IOException
	{

		long time = Calendar.getInstance().getTimeInMillis();

		dtMolecules.reset();

		setStatus("Starting processing " + dtMolecules.getRowsSize() + " memory: " + Integer.toString((int) (usedMemory() / (1024 * 1024))) + "MB");

		for (int i = 0; i < dtMolecules.getRowsSize() && dtResults.getRowsSize() != dtMolecules.getRowsSize(); i += batch)
		{
			int mol = batch != 0 && batch < (dtMolecules.getRowsSize() - i) ? batch : (dtMolecules.getRowsSize() - i);
			try
			{
				setStatus("processing " + (dtMolecules.currentRow + 1) + " out of " + dtMolecules.getRowsSize() + " molecules " +
						(mol > 1 ? "in a batch of " + mol : "") + " memory: "
						+ Integer.toString((int) (usedMemory() / (1024 * 1024))) + " MB time: " + (Calendar.getInstance().getTimeInMillis() - time) / 1000 + " sec");
				dtMolecules.setCurrentRow(i);
				calculateDescriptors(dtMolecules, receivedConfiguration, i, mol); // first attempt is to calculate all molecules in one step
			}
			catch (Exception e)
			{
				e.printStackTrace(out);
				if (batch == 1)
				{ 
					// we are working with one-by-one already: there is no need to repeat calculations
					if (dtResults.getRowsSize() != (i + 1))
						dtResults.addRow();
					dtResults.getCurrentRow().setError(e.getMessage());
					out.println("Exception " + e.getMessage());
					continue;
				}
			}

			if (batch == 1 || dtMolecules.getRowsSize() == 1)  // only one molecule was analysed
			{ 
				// we are working with one-by-one already: there is no need to repeat calculations
				if (dtResults.getRowsSize() != (i + 1)) {
					dtResults.addRow();
					dtResults.getCurrentRow().setError("calculation failed for this molecule in one by one mode");
					continue;
				}
			}

			if (dtResults.getRowsSize() != (i + mol))
			{ 
				if (getBatchSize() == 0) { // mark all calculations as errors and stop
					for (int j = 0; j < dtResults.getRowsSize() ; j++)
						dtResults.getRow(i).setError("calculation failed");

					for (int j = dtResults.getRowsSize(); j < dtMolecules.getRowsSize() ; j++)
						dtResults.addRow().setError("calculation failed");
					break;
				}

				// failed, calculate one by one
				setStatus("result " + (dtResults.getRowsSize()) + " != " + (i + mol));

				dtResults.setSize(i); // we reset the row counter and repeat calculations again
				for (int j = i; j < (i + mol); j++)
				{
					try
					{
						dtMolecules.setCurrentRow(j);
						setStatus("processing one-by-one " + (dtMolecules.currentRow + 1) + " molecule out of "
								+ dtMolecules.getRowsSize() + ". Time spent: " + (Calendar.getInstance().getTimeInMillis() - time) / 1000 + " sec");
						calculateDescriptors(dtMolecules, receivedConfiguration, j, 1);
						out.println("Intermediate result are ready for molecules:"+dtResults.getRowsSize()+" columns: "+dtResults.getColumnsSize());
						if(dtResults.getRow(j).isError())dtResults.getRow(j).setError(dtResults.getRow(j).detailedStatus);
					} catch (Exception e)
					{
						if (dtResults.getRowsSize() != (j + 1))
							dtResults.addRow();
						dtResults.getCurrentRow().setError(e.getMessage());
						out.println("Exception " + e.getMessage());
						e.printStackTrace(out);
					}
				}
			}
		}

		// If no molecules succeeded, fail completely -- not correct behavior if other molecules are cached!
		//if (!anyValidResults)
		//	throw new UserFriendlyException("Descriptor for all molecules failed. Exemplary error message: " + dtResults.getRow(0).detailedStatus);

		out.println("Result are ready for molecules: "+dtResults.getRowsSize()+" columns: "+dtResults.getColumnsSize());

		return dtResults;
	}

	/**
	 * Abstract method which should be implemented by each Descriptor Calculation Server
	 * @param dtMolecules
	 * @param configuration
	 * @param start
	 * @param batchSize
	 * @return
	 * @throws Exception
	 */

	protected abstract DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration configuration, int start, int batchSize) throws Exception;


}
