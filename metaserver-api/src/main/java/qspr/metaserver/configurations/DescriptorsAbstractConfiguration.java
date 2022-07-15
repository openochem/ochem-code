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

import javax.servlet.http.HttpServletRequest;

import com.eadmet.util.ObjectCloner;

public class DescriptorsAbstractConfiguration implements Serializable{

	private static final long serialVersionUID = 1L;
	public Integer repostSize;
	public Integer moleculeTimeout; // in minutes

	public static final Integer DEFAULT_TIMEOUT_IN_MINUTES = 10;

	DescriptorsAbstractConfiguration(){
	}

	/**
	 * will try to calculate the failed molecules using smaller subsets
	 * @return
	 */
	public boolean isLongCalculation()
	{
		return false;   
	}

	/**
	 * Provides a clean copy of DescriptorAbstractConfiguration for storage in DescriptorAbstractConfiguration
	 * @return
	 * @throws Exception
	 */

	public DescriptorsAbstractConfiguration getCopyCleanedForCache(){
		try{
			DescriptorsAbstractConfiguration c = (DescriptorsAbstractConfiguration)(ObjectCloner.deepCopy(this));
			c.cleanXML();
			return c;
		}catch(Exception e){}
		return this;
	}


	/**
	 * Cleans non-necessary elements for DescriptorAbstractConfiguration
	 */

	public void cleanXML(){
		moleculeTimeout = null;
		repostSize = null;
	}

	public boolean requires3D() {
		return false;
	}

	/**
	 * Function is only defined for descriptors which provides computation of all values for each molecule
	 * Examples are: Dragon, CDK
	 * Negative examples are: Estate, Molprint, Fragmentor, etc. which generate very large numbers of descriptors 
	 * @throws IOException 
	 *  
	 */

	public DescriptorsAbstractConfiguration configurationWithAllOnDescriptors() throws IOException{
		DescriptorsAbstractConfiguration conf = getCopyCleanedForCache();
		if(requires3D())conf.setAllOn();
		else
			conf.setAllOn2D();
		return conf;
	}

	@Override
	public String toString() {
		return "this configuration is outdated and should not be used anymore";
	}

	/**
	 * Should return not null value if operation "ON" makes sense for the particular descriptor set
	 * 1) all descriptors are be calculated (e.g., only possible for those supporting calculation of "full set of descriptors"
	 * 2) if all descriptors are ON, it is different from the default configuration (i.e, operation "ON" does make sense for it, e.g. Dragon, CDK but not Estate)
	 * @return
	 */

	public DescriptorsAbstractConfiguration setAllOn() throws IOException{
		throw new IOException("setAllOn operation is not defined for " + toString());
	}

	public Integer getAllDescriptorsNumber(){
		return null;
	}

	protected void setAllOn2D() throws IOException {
		throw new IOException("setAllOn operation is not defined for " + toString());
	}

	public void setTimeout(Integer value) {
		if(value < 1) moleculeTimeout = 1;
		else
			if(value > 99)moleculeTimeout = 99;
			else
				moleculeTimeout = value;

		if(moleculeTimeout == DEFAULT_TIMEOUT_IN_MINUTES) moleculeTimeout = null;
	}

	public int getTimeoutInMinutes() {
		return moleculeTimeout == null? DEFAULT_TIMEOUT_IN_MINUTES : moleculeTimeout;
	}

	public int getTimeoutInSeconds() {
		return getTimeoutInMinutes() * 60;
	}

	/**
	 * By default all configurations are cachables i the corresponding descriptor type is cachable
	 * @return
	 */

	public boolean isCachable(){
		return false;	
	}

	/**
	 * Indicates that order of columns should be always in the same order 
	 * This is important to retrieve model predictions.
	 * If this option is on, all values, including 0, will be stored in cache. 
	 * For most of purposes the order is not important.
	 * @return
	 */

	public boolean keepColumnsInOrder() {
		return false;
	}

	/** 
	 * To be substituted, if required, by default descriptor type (and server) name
	 * @return
	 */

	public String getDefaultTypeName(){
		return this.getClass().getCanonicalName();
	}

	DescriptorsAbstractConfiguration setConfiguration(HttpServletRequest request) throws Exception{
		return this;
	}


}