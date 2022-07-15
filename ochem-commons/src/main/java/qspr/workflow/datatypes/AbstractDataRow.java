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

package qspr.workflow.datatypes;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlTransient;

import qspr.util.ClassCompressor;
import qspr.workflow.utils.QSPRConstants;

public abstract class AbstractDataRow implements Serializable, Cloneable
{
	private static final long serialVersionUID = 2L;

	@XmlAttribute
	public String status;

	@XmlTransient
	public Map<String, Serializable> attachments;

	@XmlElement
	public String detailedStatus;

	public Map <String, Serializable> getAttachmentsCopy()
	{
		if (attachments == null)
			return null;

		Map <String, Serializable> newAttachments = new HashMap<String, Serializable>();

		for (String key : attachments.keySet())
			newAttachments.put(key, attachments.get(key));

		return newAttachments;
	}


	public AbstractDataRow setError(String message)
	{
		this.status = QSPRConstants.ERROR_STATUS;
		if (message != null)
			this.detailedStatus = message.replaceAll("[^\\u0009\\u000A\\u000D\u0020-\\uD7FF\\uE000-\\uFFFD\\u10000-\\u10FFF]+", "");
		else
			this.detailedStatus = null;
		return this;
	}

	public boolean isError()
	{
		return QSPRConstants.ERROR_STATUS.equals(this.status);
	}

	public AbstractDataRow setStatus(String status)
	{
		this.status = status;
		return this;
	}



	public AbstractDataRow addAttachment(String id, Serializable object)
	{
		if (attachments == null)
			attachments = new HashMap<String, Serializable>();
		attachments.put(id, object);

		return this;
	}

	public AbstractDataRow removeAttachment(String id)
	{
		if (attachments != null)
			attachments.remove(id);

		return this;
	}

	public Serializable getAttachment(String id)
	{
		if (attachments == null)
			return null;
		return attachments.get(id);
	}

	public AbstractDataRow getDeeperCopy() {
		return (AbstractDataRow) ClassCompressor.byteToObject(ClassCompressor.objectToByte(this));
	}

	public abstract void compact();

	/**
	 * Provide representation of the DataRow to be used, e.g., in hash calculation
	 */

	@Override
	public String toString() {
		return attachments == null? "" : " " + attachments.toString();
	}

	// TODO Check whats the difference with setWidth()
	protected abstract void pad(int size);

	public abstract void setValue(int col, Serializable value);

	public abstract Serializable getValue(int col);

	/**
	 * Number of elements in the row. Depends on type. For sparse format it is number of non zero elements
	 * @return
	 */
	public abstract int size();

	/**
	 * Add columns to the row starting with a specified position
	 * Overrides existing columns (if any)
	 * @param source
	 * @param startingColumn
	 */

	public void addColumns(AbstractDataRow source,int startingColumn)
	{
		if (source.isError())
		{
			setError(source.detailedStatus);
			return;
		}

		for (int i=0; i<source.size(); i++)
			setValue(startingColumn+i, source.getValue(i));

		if (status == null) //Merge status
		{	
			status = source.status;
			detailedStatus = source.detailedStatus;
		}
	}

	public abstract void setWidth(int columns);

	public boolean hasNAN() {
		for (int i=0; i<size(); i++)
			if(Double.isNaN((Double)getValue(i)))return true;
		return false;
	}


	public boolean sameAs(AbstractDataRow rowp) {
		return false;
	}

	public int getHashValues() {
		StringBuffer buf = new StringBuffer();
		for(int i = 0; i<size();i++)
			buf.append("_"+getValue(i));
		return buf.toString().hashCode();
	}


	public void normalise(double toValue) {
		Double sum = 0.;
		for(int i=0;i<size();i++) {
			sum += Double.valueOf(""+getValue(i)); 
		}
		if (sum == 0) return;
		sum = toValue/sum; // normalization
		for(int i=0;i<size();i++) {
			double val = Double.valueOf(""+getValue(i)); 
			if(val == 0)continue;
			setValue(i, val*sum);
		}

	}

}
