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

package com.eadmet.descriptorcache;

import java.io.StringWriter;

import javax.persistence.Transient;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.annotation.XmlRootElement;

import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorType;
import qspr.metaserver.cs.WorkflowNodeServer;

import com.eadmet.utils.OCHEMUtils;

@XmlRootElement(name = "desc-config-entry")
public class DescriptorConfigEntry
{
	public String objectID;

	public String description;
	public String md5;
	public String type;
	public String user;

	/**
	 * Indicates whether the descriptors require 2D or 3D conformation
	 */
	@Transient
	public boolean twoD = false; 

	/**
	 * Indicates that order of descriptors names is important and thus all values should be stored
	 * Namely used for storage of model predictions
	 */

	public void setUser(String user){
		this.user = user;
	}

	@Transient 
	public boolean keepAllValues = false;

	/**
	 * Calculated field. 
	 * Number of entries in the cache for this config and a particular user
	 */
	public long entriesCount;
	public long entriesSize;

	public DescriptorConfigEntry(){
	}

	public DescriptorConfigEntry(DescriptorsAbstractConfiguration conf) throws JAXBException{
		init(conf, conf.getDefaultTypeName());
	}

	public DescriptorConfigEntry(DescriptorsAbstractConfiguration conf, String descriptorType) throws JAXBException{
		init(conf, descriptorType);
	}

	public DescriptorConfigEntry(DescriptorType deskType) throws JAXBException{
		init(deskType.configuration, deskType.type);
	}

	private void init(DescriptorsAbstractConfiguration conf, String descriptorType) throws JAXBException{
		type = descriptorType;

		description = getXML(conf);
		updateMD5();
		setSpecificOption(conf);
	}


	void setSpecificOption(DescriptorsAbstractConfiguration configuration){
		twoD = configuration == null || !configuration.requires3D();
		keepAllValues = configuration != null && configuration.keepColumnsInOrder();
	}


	public String toString()
	{
		return "" + user + ":" + type + ":" + md5;
	}

	private String getXML(DescriptorsAbstractConfiguration conf) throws JAXBException
	{
		if (conf == null) return "";

		conf = conf.getCopyCleanedForCache();

		StringWriter writer = new StringWriter();
		Marshaller m = WorkflowNodeServer.jContext.createMarshaller();
		m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, new Boolean(true));
		m.marshal(conf, writer);
		return writer.toString();
	}

	public void updateMD5()
	{
		String s = type + description.replaceAll("[\\n\\s]+", "").replaceAll("<?xml[^>]*?>", "");
		md5 = OCHEMUtils.getMD5(s);
	}
}
