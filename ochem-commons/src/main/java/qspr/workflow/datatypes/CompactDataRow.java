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
import java.io.Serializable;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;

@XmlRootElement(name = "crow")
public class CompactDataRow extends AbstractDataRow {
	private static final long serialVersionUID = 2L;

	@XmlTransient
	protected TFloatArrayList value;

	@XmlElement(name = "value")
	protected byte[] getValue() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DataOutputStream dos = new DataOutputStream(baos);
		dos.writeInt(value.size());
		for (int i = 0; i < value.size(); i++)
			dos.writeFloat(value.get(i));
		dos.flush();
		return baos.toByteArray();
	}

	public float[] toArray(){
		return value.toArray();
	}

	protected void setValue(byte[] data) throws IOException, ClassNotFoundException {
		ByteArrayInputStream bais = new ByteArrayInputStream(data);
		DataInputStream dis = new DataInputStream(bais);
		int size = dis.readInt();
		value.ensureCapacity(size);
		for (int i = 0; i < size; i++)
			value.add(dis.readFloat());
		dis.close();
	}

	public CompactDataRow(int size) {
		value = new TFloatArrayList(size);
		setWidth(size);
	}

	public CompactDataRow() {
		value = new TFloatArrayList();
	}

	@Override
	public String toString() {
		return value.toString() + super.toString();
	}

	@Override
	public AbstractDataRow setError(String message) {
		value = new TFloatArrayList(0); // Empty the row
		return super.setError(message);
	}

	protected float serializebleToValue(Serializable s) {
		if (s == null)
			return 0F;

		if (s instanceof Double)
			return ((Double) s).floatValue();
		if (s instanceof Float)
			return (Float) s;
		if (s instanceof Integer)
			return (Integer) s;

		return new Float(s.toString());
	}

	@Override
	protected void pad(int size) {
		while (value.size() < size)
			value.add(0);
	}

	@Override
	public void setValue(int col, Serializable value) {
		pad(col + 1);
		this.value.set(col, serializebleToValue(value));
	}

	@Override
	public Serializable getValue(int col) {
		if (this.value.size() <= col)
			return null;
		else
			return new Double(this.value.get(col));
	}

	@Override
	public int size() {
		return value.size();
	}

	@Override
	public void compact() {
		value.trimToSize();
	}

	@Override
	public void setWidth(int columns) {
		value.ensureCapacity(columns);
		for (int i = value.size(); i < columns; i++)
			value.add(0);
	}

	public boolean sameAs(AbstractDataRow rowp) {
		if(isError() && rowp.isError()) return true;
		int size = size();
		if(size != rowp.size())return false;
		for(int i = 0; i< size; i++) {
			double a = (Double) getValue(i);
			double b = (Double) rowp.getValue(i);

			if(a != b) return false;
		}
		return true;
	}

}
