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

import java.io.IOException;
import java.sql.Timestamp;
import java.util.Calendar;

import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

@XmlRootElement
public class CacheEntry
{
	public transient boolean isNew;
	private transient Boolean sendForCalculation;

	public String objectID;
	public String user;

	protected byte[] zippedNames;
	protected byte[] zippedValues;

	/**
	 * Stores the MD5 of the SDF or the original ID of the molecule
	 */
	public String moleculeMD5;
	
	/**
	 * Stores the OCHEM ID of the molecule
	 */
	public Integer mp2;

	@XmlTransient
	public Timestamp dateCreated = new Timestamp(Calendar.getInstance().getTimeInMillis());;

	public String error;
	public int attempts;

	public final int MAXATTEMPTS = 3;

	public DescriptorConfigEntry config;

	/**
	 * Do we have to calculate this entry?
	 * We do if (a) its new or (b) it has failed less than 3 times
	 * @return
	 */
	public boolean sendForCalculation()
	{
		// Always return the same value
		if (sendForCalculation != null)
			return sendForCalculation;
		return sendForCalculation = (mp2 != null && (isNew || (attempts < MAXATTEMPTS && error != null)) || zippedNames == null);
	}

	public void doSend()
	{
		sendForCalculation = true;
	}

	public void doNotSend()
	{
		sendForCalculation = false;
	}

	public String[] getNames()
	{
		return CustomSerializer.bytesToStringArray(zippedNames);
	}

	public float[] getValues()
	{
		return CustomSerializer.bytesToFloatArray(zippedValues);
	}

	public void setNamesAndValues(String[] names, float[] values, boolean keepAllValues) throws IOException
	{
		boolean[] skip = new boolean[values.length];
		for (int i=0; i<values.length; i++)
			if (values[i] == 0 && !keepAllValues)
				skip[i] = true;
		zippedNames = CustomSerializer.stringArrayToBytes(names, skip);
		zippedValues = CustomSerializer.floatArrayToBytes(values, skip);
	}

	public void setMD5(String md5) {
		moleculeMD5 = md5;		
	}
}
