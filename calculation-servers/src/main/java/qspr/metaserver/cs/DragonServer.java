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
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import qspr.metaserver.configurations.DescriptorsAbstractDragonConfiguration;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.utils.FileUtils;

public class DragonServer extends DescriptorsAbstractExecutableServer
{
	public static final String dragonlist = "list.txt";
	public static final String dragonlog = "dragon.log";
	public static final String molFile = "molecules.sdf";
	public static final String dragoninput = "input.dps";
	public static final String dragonresult = "results.txt";

	public Integer blocks;

	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration,int start,int size) throws Exception
	{
		if (receivedConfiguration == null
				|| !(receivedConfiguration instanceof DescriptorsAbstractDragonConfiguration))
			throw new Exception(
					"Invalid configuration passed, should be instance of DragonConfiguration");
		DescriptorsAbstractDragonConfiguration configuration = (DescriptorsAbstractDragonConfiguration) receivedConfiguration;

		createDragonCfg(configuration.dragonBlocks, blocks);
		saveMolecules(dtMolecules, molFile,QSPRConstants.SDF,start,size);
		FileUtils.saveStringToFile(molFile, getAliasedFileName(dragonlist));

		String[] commands =	new String[] {getExeFile(), dragoninput, dragonlog };
		// 30 seconds per molecule on average; maximum dead-time will be 50 minutes or 
		// 5 minutes per a problematic molecule for one-by-one processing
		executeBinary(commands, dragonresult,size==1?300:size*30); 
		return readResults(dragonresult, dragonlog, size);
	}

	@Override
	int getBatchSize(){
		return 50;
	}

	void createDragonCfg(long dragonGroups,int groups) throws IOException
	{
		String data = "DRAGON script Ver 2.1\n";
		for (int i = 1; i <= groups; i++)
		{
			boolean useThisBlock = (dragonGroups % 2 == 1);
			dragonGroups /= 2;
			data += "/d GetB" + i + ((useThisBlock) ? " All" : " None")
					+ " /PCno\n";
		}
		data += "/fm list.txt -f4 -i2 -Hy\n";
		data += "/fy None\n";
		data += "/fo " + dragonresult + " -f1 -k -m -999";
		FileUtils.saveStringToFile(data, getAliasedFileName(dragoninput));
	}

	Map<String, String> getErrorsFromLog(String logFile)
	{
		Map<String, String> errors = new HashMap<String, String>();

		if(logFile == null) return errors;
		File f = new File(getAliasedFileName(logFile));
		if(!f.exists())return errors;

		try {
			BufferedReader br = getAliasedBufferedReader(logFile);
			Pattern pattern = Pattern.compile(".*Molecule (\\d+).*rejected.*",Pattern.DOTALL | Pattern.CASE_INSENSITIVE);
			String s = null;
			while ((s = br.readLine()) != null)
			{
				Matcher m = pattern.matcher(s);
				if (!m.matches())
					continue;
				try {
					String num = m.group(1).trim();
					s = br.readLine();

					if (s.contains("cause:")) //Dragon6 special case
						s = br.readLine();

					String error = s.replaceAll(".*Error: ", "");
					errors.put(num, error);
				} catch (Exception e)
				{
					//Could not properly extract error, ignore it
				}
			}
			br.close();
		} catch (Exception e)
		{
			e.printStackTrace(out);
		}
		if (errors.keySet().size() > 0)
			out.println("Managed to parse out "+errors.keySet().size()+" error messages from the log for the following molecule ids: "+errors.keySet());
		else
			out.println("No errors parsed from log file");
		return errors;
	}

	private void setErrorForLine(String dragonLine, Map<String, String> errors, AbstractDataRow currentRow)
	{
		String[] pieces = dragonLine.split("\t", 2);
		String key = pieces[0].trim();
		String error = errors.get(key);
		if (error != null)
		{
			currentRow.setError("Dragon failed: " + error);
			out.println("Set detailed error for row '"+key+"': "+error);
		}
		else
		{
			currentRow.setError("Dragon failed for unknown reasons for this molecule");
			out.println("Set generic error for row '"+key+"'");
		}
	}

	DataTable readResults(String resultFile, String logFile, int totalMolecules) throws Exception
	{
		Map<String, String> errors = getErrorsFromLog(logFile);
		BufferedReader br = getAliasedBufferedReader(resultFile);

		DataTable dtResults = getResults();

		String dragonLine = null;

		do {
			dragonLine = br.readLine();
		} while ((dragonLine != null) && (dragonLine.indexOf("No.") == -1)); // Search for line with names of columns -- it is different for different versions of Dragon

		if (dragonLine == null)
			throw new Exception("Dragon failed: reached end of file when more molecules are expected");

		String[] dragonCells = dragonLine.split("\t");

		if (dtResults.getColumnsSize() == 0)
			for (int i = 2; i < dragonCells.length; i++) // The dragon header looks like: "No. NAME Desc1 Desc2 ...."
				dtResults.addColumn(dragonCells[i]);

		for (int molecule = 0; molecule < totalMolecules; molecule++)
		{
			dtResults.addRow(); // add new Row
			dragonLine = br.readLine();

			if (dragonLine == null)
				throw new Exception("Dragon failed: reached end of file when more molecules are expected");

			if (dragonLine.contains("Error"))
			{
				setErrorForLine(dragonLine, errors, dtResults.getCurrentRow());
				continue;
			}

			dragonCells = dragonLine.split("\t");

			int normalDescriptors = 0;

			// We are fine, no error here
			for (int i = 0; i < dtResults.getColumnsSize(); i++)
			{
				if (dragonCells[i + 2].equals("NaN") || dragonCells[i + 2].equals("-999"))
					dtResults.setValue(i, 0D); // replace missing values with zero
				else
				{dtResults.setValue(i, new Double(dragonCells[i + 2]));normalDescriptors++;}
			}

			if (normalDescriptors == 0)
				setErrorForLine(dragonLine, errors, dtResults.getCurrentRow());

		}

		br.close();
		return dtResults;
	}


}
