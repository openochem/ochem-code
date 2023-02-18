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

import org.apache.commons.lang.StringUtils;

import qspr.metaserver.configurations.AbstractWekaConfiguration;
import qspr.metaserver.configurations.LabelWeighting;
import qspr.metaserver.configurations.ModelAbstractConfiguration;
import qspr.metaserver.cs.util.DescriptorsTable;
import qspr.metaserver.cs.util.LabelsTable;
import qspr.workflow.utils.QSPRConstants;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.NumericalValueStandardizer;

public abstract class WekaAbstractServer extends MachineLearningExecutableAbstractServer
{
	protected String className = "";

	private final static String cfgFileName = "cfg.xml";
	private final static String trainFileName = "train.arff";
	private final static String testFileName = "test.arff";
	private final static String modelFileName = "model.bin";
	private final static String classPrefix = "class";

	private final String HEADERS = "<!DOCTYPE options\n" + "[\n<!ELEMENT options (option)*>\n"
			+ "<!ATTLIST options type CDATA \"classifier\">\n" + "<!ATTLIST options value CDATA \"\">\n"
			+ "<!ELEMENT option (#PCDATA | options)*>\n" + "<!ATTLIST option name CDATA #REQUIRED>\n"
			+ "<!ATTLIST option type (flag | single | hyphens | quotes) \"single\">\n]\n>\n";

	private String enumerateClasses(int numOfClasses)
	{
		StringBuffer sb = new StringBuffer();
		sb.append("{");
		for (int j = 0; j < numOfClasses; j++)
		{
			if (j != 0)
				sb.append(", ");
			sb.append(WekaAbstractServer.classPrefix + j);
		}
		sb.append("}");
		return sb.toString();
	}

	protected abstract void writeMethodSpecificCfgXML(BufferedWriter out, ModelAbstractConfiguration conf) throws Exception;

	boolean useCostMatrix(AbstractWekaConfiguration conf){
		return conf.labelWeighting != null && conf.labelWeighting.useCostMatrix();
	}

	private void writeTrainingCfgXml(AbstractWekaConfiguration conf, String cfgFileName) throws Exception
	{
		BufferedWriter out = getAliasedBufferedWriter(cfgFileName);

		out.write(HEADERS);

		if (useCostMatrix(conf))
		{
			out.write("<options value=\"weka.classifiers.meta.CostSensitiveClassifier\" type=\"class\">\n");
			out.write("<option type=\"hyphens\" name=\"W\">\n");
			out.write("<options type=\"classifier\" value=\"" + className + "\">\n");
			writeMethodSpecificCfgXML(out, conf);
			out.write("</options>\n");
			out.write("</option>\n");
		}
		else
		{
			out.write("<options value=\"" + className + "\" type=\"class\">\n");
			writeMethodSpecificCfgXML(out, conf);
		}

		out.write("<option type=\"single\" name=\"p\">0</option>\n");
		out.write("<option type=\"single\" name=\"t\">" + getAliasedFileName(trainFileName) + "</option>\n");
		out.write("<option type=\"single\" name=\"d\">" + getAliasedFileName(modelFileName) + "</option>\n");
		out.write("</options>");
		out.close();
	}

	private void writeTestCfgXml(String cfgFileName) throws Exception
	{
		BufferedWriter out = getAliasedBufferedWriter(cfgFileName);
		out.write(HEADERS);
		out.write("<options value=\"" + className + "\" type=\"class\">\n");
		out.write("<option type=\"single\" name=\"p\">0</option>\n");
		out.write("<option type=\"single\" name=\"T\">" + getAliasedFileName(testFileName) + "</option>\n");
		out.write("<option type=\"single\" name=\"l\">" + getAliasedFileName(modelFileName) + "</option>");
		out.write("</options>");
		out.close();
	}

	private void writeCostMatrix(LabelWeighting weighting) throws IOException
	{
		BufferedWriter out = getAliasedBufferedWriter("wekamodel.cost");

		int numOfClasses = weighting.propertyWeights.get(0).classesWeights.size();
		out.write("% Rows\tColumns\n");
		out.write(String.format("%d\t%d\n", numOfClasses, numOfClasses));
		out.write("% Matrix elements\n");
		for (int i = 0; i < numOfClasses; i++)
			out.write(StringUtils.join(weighting.propertyWeights.get(0).classesWeights.get(i).costMatrixWeights, "\t") + "\n");
		out.close();
	}

	private void readPredictions(DataTable results, BufferedReader in) throws Exception
	{
		in.readLine();
		in.readLine();
		in.readLine();
		in.readLine();
		in.readLine();
		String line;
		while ((line = in.readLine()) != null)
		{
			String[] pieces = line.trim().split("\\s+");

			if (pieces.length < 3)
				break;

			String predicted = pieces[2];
			if (predicted.contains(WekaAbstractServer.classPrefix))
			{
				pieces = predicted.split(":");
				predicted = pieces[1].replaceAll(WekaAbstractServer.classPrefix, "");
			}
			results.addRow();
			results.setValue(new Double(predicted));
		}
	}

	private void writeHeadersARFF(DescriptorsTable dtDesc, LabelsTable dtVals,  BufferedWriter out, boolean compactDescriptors)
			throws IOException
	{

		out.write("@RELATION wekamodel\n\n");

		for (int i = 0; i < dtDesc.getRawColumnsNumber(); i++)
		{
			out.write("@ATTRIBUTE\t");
			if(compactDescriptors)out.write("\"" + i + "\"\t"); // write column Id
			else
				out.write("\"" + dtDesc.getRawColumnName(i).replaceAll("\\\\","_") + "\"\t"); // write name
			if (dtDesc.optionSize(i)!=1) // = qualitative
				out.write(enumerateClasses(dtDesc.optionSize(i)+1));
			else
				out.write("NUMERIC");
			out.write("\n");
		}

		//Weka doesn't support multilearning out of the box, pretend we have single learning for now

		out.write("@ATTRIBUTE\tclass\t");
		if (dtVals.getNumberOfClasserPerProperty(0)>1) // = classification
			out.write(enumerateClasses(dtVals.getNumberOfClasserPerProperty(0)));
		else
			out.write("NUMERIC");

		out.write("\n\n@DATA\n");
	}


	private void writeDescriptorLineARFF(DescriptorsTable dtDescriptors, int molecule, BufferedWriter out)
			throws IOException
	{
		for (int i = 0; i < dtDescriptors.getRawColumnsNumber(); i++)
		{
			double value = (Double) dtDescriptors.getRawColumnValue(molecule,i);
			int roundedValue = (int)Math.rint(value);
			if (dtDescriptors.optionSize(i)!=1)
				out.write(WekaAbstractServer.classPrefix + roundedValue + ",");
			else
				out.write(NumericalValueStandardizer.formattedFloatValuesForWEKA(value)+",");
		}
	}

	private void writeTrainingSetARFF(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, String trainFileName) throws IOException
	{
		BufferedWriter out = getAliasedBufferedWriter(trainFileName);

		writeHeadersARFF(dtDescriptors, dtExpValues, out, true);

		for(int mol=0;mol<dtExpValues.getDataSize();mol++)
		{
			writeDescriptorLineARFF(dtDescriptors, mol, out);
			if (dtExpValues.getNumberOfClasserPerProperty(0) > 1)
				out.write(WekaAbstractServer.classPrefix);
			out.write(dtExpValues.getOneValueOnlyString(mol)+"\n");
		}

		out.close();

	}

	private void writeTestSetARFF(DescriptorsTable dtDescriptors, LabelsTable dtExpValues, String testFileName, boolean compactDescriptors)
			throws IOException
	{
		BufferedWriter out = getAliasedBufferedWriter(testFileName);
		writeHeadersARFF(dtDescriptors, dtExpValues, out, compactDescriptors);

		for(int mol=0;mol<dtDescriptors.getDataSize();mol++)
		{
			writeDescriptorLineARFF(dtDescriptors, mol, out);
			out.write("?\n");
		}
		out.close();
	}

	@Override
	protected DataTable applyModelBatch(DescriptorsTable dtDescriptors,  ModelAbstractConfiguration receivedConf) throws Exception{
		AbstractWekaConfiguration wekaConf = (AbstractWekaConfiguration) receivedConf;

		LabelsTable dtExpValues=new LabelsTable(null, receivedConf);

		writeTestSetARFF(dtDescriptors, dtExpValues, testFileName, wekaConf.useCompactDescriptorNames == null ? false:wekaConf.useCompactDescriptorNames);
		writeTestCfgXml(cfgFileName);

		saveModelToFile(receivedConf, getAliasedFileName(modelFileName));

		executeBinary(getCommands(wekaConf), ExecutableServer.stdout, 0, null,null);

		DataTable dtResults = new DataTable(true);
		dtResults.addColumn(QSPRConstants.PREDICTION_RESULT_COLUMN);

		BufferedReader resultsReader = getAliasedBufferedReader(ExecutableServer.stdout);
		readPredictions(dtResults, resultsReader);
		resultsReader.close();

		return dtResults;
	}

	/**
	 * On Java we use the default java available in the path
	 * Moreover we use default 512MB 2GB is too much
	 * @return
	 * @throws IOException 
	 */

	private String[] getCommands(AbstractWekaConfiguration conf) throws IOException{

		boolean useCostMatrix = useCostMatrix(conf);
		if (!conf.isTrainingConfiguration())
			useCostMatrix = false;

		String memory =  "-Xmx"+ getMemoryForExecutable()+"M";

		if (!useCostMatrix)
			return new String[] {javaHome+"/bin/java", "-cp", javaClassPath, memory, className, "-xml", cfgFileName };
		else
		{
			writeCostMatrix(conf.labelWeighting);
			return new String[] {javaHome+"/bin/java", "-cp", javaClassPath, memory, "weka.classifiers.meta.CostSensitiveClassifier", "-xml", cfgFileName};
		}
	}

	@Override
	protected DataTable trainModel(DescriptorsTable dtDescriptors, LabelsTable dtExpValues,
			ModelAbstractConfiguration receivedConf) throws Exception
	{		

		AbstractWekaConfiguration wekaConf = (AbstractWekaConfiguration) receivedConf;
		wekaConf.useCompactDescriptorNames = true;

		checkMethodSpecificPrerequisites(wekaConf);

		writeTrainingSetARFF(dtDescriptors, dtExpValues, trainFileName);
		writeTrainingCfgXml(wekaConf,cfgFileName);

		executeBinary(getCommands(wekaConf), ExecutableServer.stdout, 0, null,null);

		receivedConf.storeModel(getAliasedFile(modelFileName));

		return applyModel(dtDescriptors, wekaConf);
	}

	/**
	 *  Check which type of tasks are supported
	 * @param dtDescriptors
	 * @param dtExpValues
	 * @param receivedConf
	 * @throws Exception
	 */
	protected void checkMethodSpecificPrerequisites(ModelAbstractConfiguration receivedConf) throws Exception
	{
		if (receivedConf.areMultiLearningData())
			throw new UserFriendlyException(supportedTaskType+ " doesn't support multilearning. Please, provide a dataset with single property only.");

		if (!receivedConf.areClassificationData())
			throw new UserFriendlyException(supportedTaskType + " doesn't support regression tasks. Please, provide a classification task.");

	}

}
