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

package qspr.entities;

import java.util.HashMap;
import java.util.Map;

import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import qspr.Globals;
import qspr.interfaces.Descriptable;
import qspr.metaserver.configurations.StandartizationOptions;

@XmlRootElement(name = "model-attachment")
public class ModelAttachment implements Descriptable
{
	@XmlElement(name="configuration")
	public Object configuration;
		
	public Map<Long, Long> optionsMapping = new HashMap<Long, Long>(); // option id -> qualitative value
	
	public ModelProtocol protocol = new ModelProtocol();
	public StandartizationOptions standartization = new StandartizationOptions();
	public DataHandlingOptions datahandling = new DataHandlingOptions();
		
	//@XmlElement(name = "ad-configuration")
	//public ADConfiguration adConfiguration;
	
	public void printXml() throws JAXBException
	{
		Marshaller marshaller = Globals.jaxbContext.createMarshaller();
		marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT,
				   new Boolean(true));
		marshaller.marshal(this, System.out);
	}
	
	/**
	 * Show hierarchical configuration of this model as a string / Midnighter
	 */
	
	public String toString()
	{
		return "" + (configuration != null ? configuration + "\n" : "") + protocol;
	}
	
	public Map<Long, Long> getOptionsMapping()
	{
		Map<Long, Long> options = new HashMap<Long, Long>();
		options.putAll(optionsMapping);
		return options;
	}
	
	public Map<String, Object> getParameters() 
	{
		Map<String, Object> parameters = new HashMap<String, Object>();
		if (configuration instanceof Descriptable)
			parameters.putAll(((Descriptable)configuration).getParameters());
		
		parameters.put("Validation", (protocol.validationConfiguration == null ? "No validation" : ""+protocol.validationConfiguration.ensembleSize+" CV"));
		return parameters;
	}
	
	public PropertyOption getOptionFromPrediction(Double predictionValue, Property property)
	{
		predictionValue = Double.valueOf(Math.round(predictionValue));
		Long closestLabel = 999999999L; // in case we have a double value
		PropertyOption result = null;
		for (Long optionId : optionsMapping.keySet()) 
		{
			PropertyOption po = (PropertyOption) Globals.session().get(PropertyOption.class, optionId);
			if (po.property.equals(property))
				if(predictionValue.equals(optionsMapping.get(optionId).doubleValue()))
					return po;
				else if (Math.abs(predictionValue - closestLabel) > Math.abs(predictionValue - optionsMapping.get(optionId)))
				{
					closestLabel = optionsMapping.get(optionId);
					result = po;
				}
		}
		
		return result;
	}	
	
}
