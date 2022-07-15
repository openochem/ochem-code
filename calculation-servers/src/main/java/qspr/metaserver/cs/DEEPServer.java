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

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.Writer;

import com.eadmet.utils.FileUtils;
import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.OSType;

import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.configurations.MultiLearningAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;


public class DEEPServer extends MultiLearningAbstractServer
{

	private static final String DATAFILE = "boston";
	private static final String MODEL = "model";

	private static final String ANN="ANN.dbd";
	
	int outputs =0;
	
	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws Exception
	{
		//create data file for Training
		outputs= receivedConf.optionsNumber == null?1:receivedConf.OutputValues();
		saveAggregatedData(DATAFILE, dtDescriptors, dtExpValues, receivedConf);

		if(outputs>1) {
			setStatus("Analysing " + receivedConf.OutputValues() + " " + dtExpValues.getColumnSize());
			FileUtils.saveStringToFile("  "+outputs, getAliasedFileName(ANN));
		}
		
		String[] commands = new String[] { getExeFile(), DATAFILE, outputs>1?"--SCALE_PAT_MUL":"--SCALE_DEEP"};
		executeBinary(commands, DATAFILE+".txt", 0);
		commands[1] = DATAFILE+".txt"; commands[2] = "--SPLITT_BIG_A";
		executeBinary(commands, DATAFILE+".pat", 0);
		commands[1] = DATAFILE+".pat"; commands[2] = outputs>1?"--NET5":"--NET4";
		executeBinary(commands, DATAFILE+".dbd", 0);
		commands[1] = DATAFILE+".pat"; commands[2] = "--BP";

		executeBinary(commands, DATAFILE+".wgt", 0);

		commands[0] = "tar";
		commands[1] = "-cf";
		commands[2] = MODEL;
		commands = OCHEMUtils.append(commands, DATAFILE+".wgt");
		commands = OCHEMUtils.append(commands, DATAFILE+".dbd");
		if(outputs==1)commands = OCHEMUtils.append(commands, "la_scala.txt");
		else
			commands = OCHEMUtils.append(commands, ANN);
		commands = OCHEMUtils.append(commands, "statts.txt");

		executeBinary(commands, MODEL, 0);

		receivedConf.storeModel(getAliasedFile(MODEL));

		setStatus("DEEP has finished; applying the model");
		return applyModel(dtDescriptors, receivedConf);
	}


	public DEEPServer(){
		supportedTaskType = QSPRConstants.DEEP;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
		DELIMITER = "\\s+";
	}

	@Override
	protected DataTable applyModelMulti(DescriptorsTable dtDescriptors, MultiLearningAbstractConfiguration receivedConf)
			throws Exception {
		BufferedWriter writer = getAliasedBufferedWriter(DATAFILE);
		outputs= receivedConf.optionsNumber == null?1:receivedConf.OutputValues();
		saveAggregatedData(DATAFILE, dtDescriptors, null, receivedConf);
		writer.close();

		saveModelToFile(receivedConf, getAliasedFileName(MODEL));
		String untar[]=new String[] { "tar","-xf",MODEL};
		executeBinary(untar, "statts.txt", 0);
		String[] commands = new String[] { getExeFile(), DATAFILE, outputs>1?"--SCALE_TES_MUL":"--SCALE_DEEP_TEST"};
		executeBinary(commands, DATAFILE+".txt", 0);

		commands[1] = DATAFILE+".txt"; commands[2] = "--BP_TEST";
		executeBinary(commands, "resultss.ttt", 0);
		commands[1] = "resultss.ttt"; commands[2] = outputs>1?"--DESCALE_MUL":"--DESCALE_DEEP";
		executeBinary(commands, "results_H.ttt", 0);
		return readResultValues("results_H.ttt",receivedConf);
	}

	@Override
	protected void saveOneRow(int mol, DescriptorsTable dtDescriptors, String[] values, Writer writer)
			throws IOException {

		String 	VALUESSEPARATOR =",  ";

		String []desc = dtDescriptors.getScaledValuesString(mol);

		writer.append("  "+mol+VALUESSEPARATOR);

		for(int i=0;i<desc.length;i++)
			writer.append(desc[i]+VALUESSEPARATOR);

		if(values != null)
			for(int i=0; i <outputs; i++) {
				if(values[i] == null || values[i].equals("null"))values[i] = "-999";
				writer.append(values[i]+ ((outputs == i+1)?"":VALUESSEPARATOR));
			}else
				for(int i=0; i <outputs; i++) 
					writer.append(0+ ((outputs == i+1)?"":VALUESSEPARATOR));

		writer.append(OSType.endLine());	
	}
		
}
