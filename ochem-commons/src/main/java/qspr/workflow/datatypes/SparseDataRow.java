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

import gnu.trove.list.array.TFloatArrayList;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.BitSet;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement(name = "srow")
public class SparseDataRow extends CompactDataRow {

	private static final long serialVersionUID = 1L;

	BitSet usedBits;

	int currentColumn=Integer.MAX_VALUE,currentPosition=0; // added to speed up work with Sparse format

	@XmlElement(name="value")
	protected byte[] getValue() throws IOException
	{
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DataOutputStream dos = new DataOutputStream(baos);
		ObjectOutputStream oos=new ObjectOutputStream(dos); 
		dos.writeInt(value.size());
		for (int i=0; i<value.size(); i++)
			dos.writeFloat(value.get(i));
		dos.flush();
		oos.writeObject(usedBits);
		oos.flush();
		oos.close();
		dos.close();
		return baos.toByteArray();
	}

	@Override
	public float[] toArray(){
		float val[]=value.toArray();
		float values[]=new float[usedBits.length()];

		for (int j=0,i = usedBits.nextSetBit(0); i >= 0; i = usedBits.nextSetBit(i+1))
			values[i]=val[j++];

		return values;
	}

	@Override
	protected void setValue(byte[] data) throws IOException, ClassNotFoundException
	{
		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		DataInputStream dis = new DataInputStream(bais);
		ObjectInputStream ois = new ObjectInputStream(dis);
		int size = dis.readInt();
		value.ensureCapacity(size);
		for (int i=0; i<size; i++)
			value.add(dis.readFloat());
		usedBits=(BitSet)ois.readObject();
		ois.close();
		dis.close();
	}	

	@Override
	public Serializable getValue(int col)
	{
		if(!usedBits.get(col))return new Double(0.); // We always return 0 for elements that do not exist
		return new Double(value.get(getPosition(col)));
	}

	/**
	 * The most critical part for performance
	 * It is optimized to get maximum speed
	 * @param col
	 * @return
	 */

	private int getPosition(int col)
	{

		if(col < currentColumn)
			currentColumn = currentPosition = 0;

		while (currentColumn < col)
		{
			if (usedBits.get(currentColumn))
				currentPosition++;
			currentColumn++;
		}
		//currentPosition = currentPosition;
		return currentPosition;
	}


	/**
	 *  This operation does not make sense for Sparse format -- we cannot allocate size in advance
	 */
	@Override
	protected void pad(int size){}

	/**
	 *  This operation does not make sense for Sparse format -- we cannot allocate size in advance
	 */

	@Override
	public void setWidth(int columns) {}


	@Override
	public void setValue(int col, Serializable newValue)
	{
		float v=serializebleToValue(newValue);
		if(usedBits.get(col)){  // we have previous non-zero value 
			int position=getPosition(col);
			if(v!=0)value.set(position,v);
			else{
				value.removeAt(position);
				usedBits.clear(col);
			}
		}else  // previous value was 0
			if(v!=0){
				usedBits.set(col);
				int position=getPosition(col);
				value.insert(position,v);
			}
	}

	public SparseDataRow()
	{
		value = new TFloatArrayList();
		usedBits=new BitSet();
	}

	@Override
	public String toString()
	{
		int size=usedBits.length();
		TFloatArrayList val=new TFloatArrayList(size);
		for(int i=0;i<size;i++)val.add(serializebleToValue(getValue(i)));
		return "data: "+val.toString()+"\nbits: "+usedBits.toString() + "\nvalues: "+super.toString()+ (isError() ? detailedStatus : "");
	}

	@Override
	public AbstractDataRow setError(String message)
	{
		value = new TFloatArrayList(0); //Empty the row
		usedBits = new BitSet(0); // Empty the bits
		currentColumn = Integer.MAX_VALUE;
		currentPosition = 0; 
		return super.setError(message);
	}

	/**
	 * In case of Sparse format we cannot know
	 * the length of a row. Thus, the length of a row
	 * is determined by the last non-zero element in the row
	 */

	@Override
	public int size()
	{
		return usedBits.length();
	}



}
