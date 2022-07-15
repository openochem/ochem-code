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

package qspr.metaserver.configurations;

import java.io.IOException;
import java.io.Serializable;

import javax.xml.bind.annotation.XmlTransient;

public  class DescriptorType implements Serializable
{

	public Boolean replaced;

	private static final long serialVersionUID = 1L;
	public String type;
	public DescriptorsAbstractConfiguration configuration;
	public String title;

	/**
	 * Should we write the results to the cache?
	 * If null, use defaults for this descriptor type
	 */
	@XmlTransient
	public Boolean writeToCache; 

	/**
	 * Should we skip the cache altogether (even reading from it)? 
	 */
	@XmlTransient
	public Boolean skipCache;

	public Boolean markUncachedAsErrors;

	public DescriptorType()
	{

	}


	public DescriptorType(DescriptorsAbstractConfiguration configuration)
	{
		this.type = configuration.getDefaultTypeName();
		this.configuration = configuration;
	}

	public DescriptorType(String taskType, DescriptorsAbstractConfiguration configuration)
	{
		this.type = taskType;
		this.configuration = configuration;
	}


	public DescriptorType(DescriptorType descr) {
		type = descr.type;
		configuration = descr.configuration;
		title = descr.title;
		writeToCache = descr.writeToCache;
		skipCache = descr.skipCache;
		markUncachedAsErrors = descr.markUncachedAsErrors;
	}

	public String toString()
	{
		if (configuration != null && configuration.toString() != "")
			return type + " ("+configuration.toString()+")";

		return type;
	}

	public boolean equals(Object obj)
	{
		return type.equals(((DescriptorType) obj).type);
	}

	/**
	 * Skip the cache for this descriptor type altogether (neither read from nor write to the cache)
	 * Only used so far for standalone calculations and Structural Alerts
	 * @return
	 */
	public DescriptorType skipCache()
	{
		skipCache = true;
		return this;
	}

	/**
	 * Skip the cache for this descriptor type altogether (neither read from nor write to the cache)
	 * @return
	 */
	public boolean isSkipCache()
	{
		return skipCache != null && skipCache;
	}

	public boolean requires3D()
	{
		return configuration != null && configuration.requires3D();
	}

	public boolean supportSplicing() {
		if(configuration == null)return false;
		try{
			configuration.configurationWithAllOnDescriptors();
		}catch(Exception e){
			return false;
		}
		return true;
	}

	public DescriptorType configurationWithAllOnDescriptors() throws IOException{
		DescriptorType c = new DescriptorType(this);
		c.configuration = c.configuration.configurationWithAllOnDescriptors();
		return c;

	}

	public boolean similar(DescriptorType a){
		return type.equals(a.type) && toString().equals(a.toString());
	}

}
