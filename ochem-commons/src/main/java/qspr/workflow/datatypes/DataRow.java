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

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlElements;
import javax.xml.bind.annotation.XmlRootElement;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

@XmlRootElement(name = "row")
public class DataRow extends AbstractDataRow
{
	private static transient final Logger logger = LogManager.getLogger(DataRow.class);

	private static final long serialVersionUID = 2L;

	@XmlElements({
		@XmlElement(name="double", type=Double.class),
		@XmlElement(name="float", type=Float.class),
		@XmlElement(name="long", type=Long.class),
		@XmlElement(name="int", type=Integer.class),
		@XmlElement(name="string", type=String.class),
		@XmlElement(name="null", type=EmptyCell.class),
		@XmlElement(name="bytearray", type = byte[].class),
		@XmlElement

	})
	//@XmlElement
	private List<Object> value;

	@XmlAttribute
	public boolean compressStrings = false;

	private Serializable preprocessValueSet(Serializable value)
	{
		if (compressStrings && value instanceof String && ((String)value).length() > 1000)
		{
			try
			{
				ByteArrayOutputStream os = new ByteArrayOutputStream();
				GZIPOutputStream gos = new GZIPOutputStream(os);
				ObjectOutputStream oos = new ObjectOutputStream(gos);
				oos.writeObject(value);
				oos.flush();
				oos.close();
				byte[] data = os.toByteArray();
				return data;
			} catch (IOException e)
			{
				logger.info("WARNING: Error in experimental String Compression routine, continuing with no compression");
				e.printStackTrace();
				return value;
			}
		}
		else	
			return value;
	}

	private Serializable preprocessValueGet(Serializable value)
	{
		if (value instanceof byte[])
		{
			try
			{
				ByteArrayInputStream is = new ByteArrayInputStream((byte[])value);
				GZIPInputStream gis = new GZIPInputStream(is);
				ObjectInputStream ois = new ObjectInputStream(gis);
				Serializable result = (Serializable)ois.readObject();
				ois.close();
				return result;
			} catch (Exception e)
			{
				//				logger.info("WARNING: Error in experimental String Decompression routine, continuing with no decompression");
				//				e.printStackTrace();
				return value;
			}
		}
		else	
			return value;		
	}

	public DataRow()
	{
		value = new ArrayList<Object>();
	}

	public DataRow(int size)
	{
		this();
		setWidth(size);
	}


	@Override
	public String toString()
	{
		return value.toString() +  super.toString();
	}

	@Override
	protected void pad(int size)
	{
		while (value.size() < size)
			value.add(new EmptyCell());
	}

	@Override
	public void setValue(int col, Serializable value)
	{
		pad(col+1);
		this.value.set(col, preprocessValueSet(value));
	}

	public void addValue(Serializable value)
	{
		this.value.add(preprocessValueSet(value));
	}	

	@Override
	public Serializable getValue(int col)
	{
		Serializable value;
		if (this.value.size() <= col)
			return null; 
		else
			//Changed to force type conversion 03.11.08
			value = (Serializable)this.value.get(col);

		if (value instanceof EmptyCell)
			value = null;

		return preprocessValueGet(value);
	}

	@Override
	public int size()
	{
		return value.size();
	}


	@Override
	public void compact()
	{
		// TODO Auto-generated method stub
		// This type of row is incompactible, I guess
	}

	@Override
	public void setWidth(int columns) 
	{
		for (int i=value.size(); i<columns; i++)
			value.add(null);
	}
}
