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
import java.io.Serializable;

import qspr.metaserver.util.ExecutableRunner;
import qspr.workflow.datatypes.DataTable;

/**
 * Executable server 
 * @author itetko
 *
 */

public abstract class ExecutableServer extends SimpleAbstractServer
{
	public final static String stdout = "stdout";
	public final static String datain = "datain";
	public final static String dataout = "dataout";

	private ExecutableRunner server;

	protected void init() throws IOException, InterruptedException{
		server = new ExecutableRunner(this);
		server.init();
	}

	protected void cleanup() throws IOException, InterruptedException{
		server.cleanup();
		server=null;
	}

	protected String getAliasedFileName(String filename) {
		return server.getAliasedFileName(filename);
	}

	protected BufferedWriter getAliasedBufferedWriter(String filename) throws IOException {
		return server.getAliasedBufferedWriter(filename);
	}

	protected BufferedReader getAliasedBufferedReader(String filename) throws IOException {
		return server.getAliasedBufferedReader(filename);
	}

	protected void executeBinary(String[] commands, String outputFile) throws IOException, InterruptedException {
		server.executeBinary(commands, outputFile);

	}

	protected void executeBinary(String[] commands, String outputFile, int timeout) throws IOException, InterruptedException {
		server.executeBinary(commands, outputFile,timeout);
	}

	protected String getExeFile() {
		return  server.getExeFile();
	}

	protected void saveMolecules(DataTable dtMolecules, String filename, String datatype, int start,int batchsize) throws IOException {
		server.saveMolecules(dtMolecules,filename,datatype,start,batchsize);		
	}


	protected void saveMolecules(DataTable dtMolecules, String filename, String datatype) throws IOException {
		server.saveMolecules(dtMolecules,filename,datatype);		
	}


	@Override
	DataTable calculateGeneralTask(DataTable dtInputData, Serializable configuration) throws Exception {
		init();
		DataTable dtResults = calculateExecutableTask(dtInputData, configuration);
		cleanup();
		return dtResults;
	}


	/**
	 *  Provides a most typical interface to calculate tasks
	 *  Usually  some data are coming and results of transformation are provided
	 * @param dtInpout
	 * @param configuration
	 * @return
	 * @throws Exception
	 */

	abstract protected DataTable calculateExecutableTask(DataTable dtInpout, Serializable configuration) throws Exception;
}
