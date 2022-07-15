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
import java.io.File;
import java.io.IOException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import qspr.metaserver.configurations.KNNConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.FileUtils;

public class KNNServer extends IgorsAbstractServer
{
	private static transient final Logger logger = LogManager.getLogger(KNNServer.class);

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws Exception
	{
		KNNConfiguration configuration = (KNNConfiguration) receivedConf;

		//create cfg file for KNN
		createCfgFile(configuration);

		//create data file for KNN Training
		BufferedWriter writer = getAliasedBufferedWriter(DATAFILE);
		writeSimpleIgorFormatData(dtDescriptors, dtExpValues, writer);
		writer.close();
		//once the cfg and data file is created then execute KNN
		executeBinary();

		setStatus("KNN has finished");
		setStatus("Reading results. Please wait");
		DataTable dtResult = readResult("ochem");
		if(dtResult.getRowsSize() != dtDescriptors.getDataSize())
			throw new Exception("Number of rows from results table did not matched with expected rows. " +
					"Number of rows in results table are :"+dtResult.getRowsSize()+
					"Number of expected rows are :"+dtDescriptors.getDataSize());

		//compress cfg, data and model and encode it for later Applier use
		setParamsToSave(configuration);
		return dtResult;
	}

	@Override
	protected DataTable applyModelBatch(DescriptorsTable dtDescriptors,  ModelAbstractConfiguration receivedConf) throws Exception
	{
		// for application, we need only results for the test set

		KNNConfiguration configuration= (KNNConfiguration)receivedConf;		
		dtDescriptors = new DescriptorsTable(configuration.getTrainingSetDescriptors().addRowsFrom(dtDescriptors.getRawData()),receivedConf,0);
		LabelsTable dtExpValues = new LabelsTable(configuration.getTrainingSetValues(), receivedConf);

		DataTable dtResult = trainModel(dtDescriptors, dtExpValues, receivedConf);
		return dtResult.getSlice(dtExpValues.getDataSize(), dtDescriptors.getDataSize());
	}

	void setParamsToSave(KNNConfiguration configuration) throws IOException, ParserConfigurationException, SAXException{
		File file = new File(getAliasedFileName(MODEL));
		if(!file.exists())throw new IOException("The "+MODEL+" file cannot be open");
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		Document doc = db.parse(file);
		doc.getDocumentElement().normalize();
		NodeList knn = doc.getElementsByTagName("KNN");
		if(knn.getLength()!=1)throw new IOException("The KNN was not specified in the model file");
		out.println("selected: "+knn.item(0).getChildNodes().item(0).getNodeValue());
		configuration.knn=(int)Long.parseLong(knn.item(0).getChildNodes().item(0).getNodeValue());
	}

	private void createCfgFile(KNNConfiguration configuration) throws IOException{
		logger.info("Writing cfg file");
		String parameters = "FILE=data\nNONZERO=1\nCORRELATION=1\nNORM=1\nNAMES=1\nMODELS=1\n";

		parameters += "\nDISTANCE="+configuration.distance;
		parameters += "\nMAXKNN="+configuration.maxKnn;
		parameters += "\nKNN="+configuration.knn;

		FileUtils.saveStringToFile(parameters+"\nSTOP\nEND", getAliasedFileName(CFG));
	}

	public KNNServer(){
		supportedTaskType = QSPRConstants.KNN;
		setInputFlowGroup(0);
		setInputFlowGroup(1);
		setOutputFlowGroup(0);
	}
}
