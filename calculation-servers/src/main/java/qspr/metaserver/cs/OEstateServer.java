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

import java.io.IOException;
import java.io.StringReader;
import java.util.List;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlValue;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsOEstateConfiguration;
import qspr.workflow.utils.QSPRConstants
;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.utils.FileUtils;

public class OEstateServer extends DescriptorsAbstractExecutableServer
{
	private static final String resultFile="indices.xml";
	private static final String cfgFile="cfg_i";
	private static final String dataFile="data.sdf";

	@Override
	int getBatchSize()
	{
		return 100;
	}
	
	@Override
	public DataTable calculateDescriptors(DataTable dtMolecules, DescriptorsAbstractConfiguration receivedConfiguration,int start,int batchsize) throws Exception
	{
		if (receivedConfiguration == null || !(receivedConfiguration instanceof DescriptorsOEstateConfiguration))
			throw new Exception("Invalid configuration passed, should be instance of OEstateConfiguration");
		DescriptorsOEstateConfiguration configuration = (DescriptorsOEstateConfiguration) receivedConfiguration;

		//create cfg file
		createCfg(configuration);
		saveMolecules(dtMolecules, dataFile,QSPRConstants.SDF,start,batchsize);	
		//create empty result table

		String[] commands = new String[]{getExeFile(), cfgFile, dataFile};
		executeBinary(commands,resultFile, batchsize==1?30:batchsize);

		//read result
		return readResult();
	}


	private void createCfg(DescriptorsOEstateConfiguration configuration) throws IOException
	{
		String data = "FILE=data.sdf\nCOUNTS="+configuration.count+"\nEXTENDED=1\nBOND="+configuration.bond+"\nXML=1\nSTANDARD=255"+
				(configuration.info == null?"":configuration.info)+			
				"\nSTOP\nEND";
		FileUtils.saveStringToFile(data, getAliasedFileName(cfgFile));
	}

	private DataTable readResult() throws Exception
	{
		String indices = FileUtils.getFileAsString(getAliasedFileName(resultFile));
		indices.replaceAll("NA>", "NA&gt;");
		StringReader sr = new StringReader(indices);
		JAXBContext context = JAXBContext.newInstance(new Class[]{EstateMolecules.class});
		Unmarshaller unmarshaller = context.createUnmarshaller();
		EstateMolecules emolecules = (EstateMolecules)unmarshaller.unmarshal(sr);

		DataTable dtResults = getResults();

		if(emolecules.molecules != null && emolecules.molecules.size() != 0)
			for(EstateMolecule emolecule :emolecules.molecules)
			{
				dtResults.addRow();
				if (emolecule.error != null)
					dtResults.getCurrentRow().setError(emolecule.error);
				else if (emolecule.descriptors == null || emolecule.descriptors.size() == 0)
					dtResults.getCurrentRow().setError("No descriptors returned by estate program");
				else
					for (EstateDescriptor edescriptor : emolecule.descriptors) 
						dtResults.setValue(edescriptor.name, edescriptor.value);
			}
		return dtResults;
	}	

	public OEstateServer()
	{
		supportedTaskType = "OEstate";
		setInputFlowGroup(0);
		setOutputFlowGroup(0);
	}
}

@XmlRootElement(name="molecules")
class EstateMolecules
{
	@XmlElement(name = "molecule")
	List<EstateMolecule> molecules;
}

class EstateMolecule
{
	@XmlAttribute
	String number;

	@XmlAttribute
	String error;

	@XmlElement(name = "descriptor")
	List<EstateDescriptor> descriptors;
}

class EstateDescriptor
{
	@XmlAttribute
	String name;

	@XmlValue
	Double value;
}
