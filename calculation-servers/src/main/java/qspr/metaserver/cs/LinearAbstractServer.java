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

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import qspr.metaserver.configurations.LinearConfiguration;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.metaserver.util.DataScaling;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;

public abstract class LinearAbstractServer extends IgorsAbstractServer {

	private static transient final Logger logger = LogManager.getLogger(LinearAbstractServer.class);

	@Override
	protected DataTable applyModelBatch(DescriptorsTable dtDescriptors, ModelAbstractConfiguration receivedConf) throws Exception 
	{
		return applyLinearModel(dtDescriptors, (LinearConfiguration) receivedConf);
	}

	/**
	 * Train model and performs normalisation of coefficients for equation 
	 */
	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws Exception 
	{
		LinearConfiguration c=(LinearConfiguration) receivedConf;
		if(c.limitRange != null && c.limitRange)limitRange(dtExpValues,c);
		return trainModelOrUseUploaded(dtDescriptors, dtExpValues, receivedConf, null);
	}

	private void limitRange(LabelsTable dtExpValues, LinearConfiguration c) throws IOException {
		c.minVal = (double)dtExpValues.getMin()[0];
		c.maxVal = (double)dtExpValues.getMax()[0];
	}

	@SuppressWarnings("unchecked")
	@Override
	protected DataTable processUploadedModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf) throws Exception 
	{
		Map<Integer, Double> coefficient = new HashMap<Integer, Double>();
		Map<String, Double> uploadedModel = (Map<String, Double>) receivedConf.uploadedModelData.get();
		for (String key : uploadedModel.keySet()) 
		{
			Double value = uploadedModel.get(key);
			if (key.equalsIgnoreCase("BIAS"))
				coefficient.put(0, value);
			else if (dtDescriptors.columnNames.contains(key))
				coefficient.put(dtDescriptors.columnNames.indexOf(key) + 1, value);
		}

		if (coefficient.get(0) == null) // If Bias is 0 or not present, we may have unwanted null here
			coefficient.put(0, 0D);

		return trainModelOrUseUploaded(dtDescriptors, dtExpValues, receivedConf, coefficient);
	}
	protected DataTable trainModelOrUseUploaded(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, ModelAbstractConfiguration receivedConf, Map<Integer, Double> uploadedCoefficients) throws Exception 
	{
		LinearConfiguration configuration = (LinearConfiguration) receivedConf;

		Map<Integer, Double> coefficient = (uploadedCoefficients == null) ? calculateCoefficients(dtDescriptors, dtExpValues, configuration) : uploadedCoefficients;

		configuration.coefficient = coefficient;

		out.println(configuration.theEquation(dtDescriptors.getRawData().getColumns(), false));

		DataScaling scale = dtDescriptors.getScaling();
		if (scale != null) 
		{
			double bias = coefficient.get(0);
			for (Integer i : coefficient.keySet()) 
			{
				if (i == 0)
					continue;
				bias += coefficient.get(i) * scale.getBias(i - 1);
				coefficient.put(i, coefficient.get(i) * scale.getSlope(i - 1));
			}
			coefficient.put(0, bias);
			dtDescriptors.setScaling(null); // forbid other normalisation
			configuration.scaleX = null;
		}

		scale = dtExpValues.getScaling();
		if (scale != null) 
		{
			double slope = scale.getInverseSlope(0);
			for (Integer i : coefficient.keySet())
				coefficient.put(i, coefficient.get(i) * slope);
			coefficient.put(0, coefficient.get(0) + scale.getInverseBias(0));

			dtExpValues.setScaling(null); // forbid other normalisation
			configuration.scaleY = null;
		}

		configuration.coefficient = coefficient;
		addNormaliseCoefficients(dtDescriptors, configuration, dtExpValues.getDataSize());

		// out.println(configuration.theEquation(dtDescriptors.getRawData().columns, false));
		out.println(configuration.theEquation(dtDescriptors.getRawData().getColumns(), true));

		DataTable res = applyModel(dtDescriptors, configuration); // using the actual values; the values will be renormalized anyway
		return res;
	}



	protected abstract Map<Integer, Double> calculateCoefficients(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, LinearConfiguration configuration) throws Exception;

	/*
	 * Get coefficients from model file
	 */
	protected Map<Integer, Double> getCoefficients(String filename, boolean fsmlra) throws IOException, ParserConfigurationException, SAXException 
	{
		Map<Integer, Double> factors = new HashMap<Integer, Double>();
		File file = new File(getAliasedFileName(filename));
		DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
		DocumentBuilder db = dbf.newDocumentBuilder();
		Document doc = db.parse(file);
		doc.getDocumentElement().normalize();
		NodeList variableLst = doc.getElementsByTagName("VARIABLE");
		NodeList factorLst = doc.getElementsByTagName("TERM");

		// get constant
		factors.put(
				0,
				Double.parseDouble(factorLst.item(0).getAttributes()
						.getNamedItem(QSPRConstants.VALUE).getNodeValue()));
		// get all the factors with their respective variables
		for (int i = 1; i < factorLst.getLength(); i++) {
			int key = Integer.parseInt(variableLst.item(i - 1).getAttributes()
					.getNamedItem("NO").getNodeValue())
					+ (fsmlra ? 1 : 0);
			double value = Double.parseDouble(factorLst.item(i).getAttributes()
					.getNamedItem(QSPRConstants.VALUE).getNodeValue());
			logger.info("variable " + key + " factor " + value);
			factors.put(key, value);
		}
		return factors;
	}

	protected DataTable applyLinearModel(DescriptorsTable dtDescriptors, LinearConfiguration configuration) throws Exception 
	{
		DataTable resultTable = new DataTable(true);
		resultTable.addColumn(QSPRConstants.PREDICTION_RESULT_COLUMN);
		int allMolecules = dtDescriptors.getDataSize();

		for (int m = 0; m < allMolecules; m++) {

			double[] v = getDescriptorValues(dtDescriptors, m);
			double res = configuration.coefficient.get(0); // A constant, should
			// be always present
			for (int i = 0; i < v.length; i++) {
				if (!configuration.coefficient.containsKey(i + 1))
					continue;
				res += configuration.coefficient.get(i + 1) * v[i];
			}
			resultTable.addRow();
			if(Double.isNaN(res))resultTable.getCurrentRow().setError("Value is NaN");
			else
			{
				if(configuration.limitRange != null && configuration.limitRange){
					if(res<configuration.minVal) res = configuration.minVal;
					if(res>configuration.maxVal) res = configuration.maxVal;
				}
				if(configuration.optionsNumber != null && configuration.optionsNumber.get(0) == 2) res = correct(res); // required to avoid non-sense values for linear methods applied to classification tasks
				resultTable.setValue(res);
			}
		}

		return resultTable;
	}

}
