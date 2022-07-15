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

import java.io.Serializable;

import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlValue;

import qspr.util.ClassCompressor;

@XmlRootElement(name = "compressed-object")
public class CompressedObject<T> implements Serializable
{
	private static final long serialVersionUID = -4808242473256224554L;

	@XmlTransient
	private byte[] compressedObject;

	private transient T cachedObject;

	@XmlValue
	protected byte[] getCompressedObject() {
		if (compressedObject == null && cachedObject != null) {
			set(cachedObject);
		}

		return compressedObject;
	}

	protected void setCompressedObject(byte[] object) {		
		this.compressedObject = object;
	}

	public void set(T obj)
	{
		cachedObject = obj;
		compressedObject = ClassCompressor.objectToByte((Serializable) obj);
	}

	@SuppressWarnings("unchecked")
	public T get()
	{
		if (cachedObject != null)
		{
			return cachedObject;
		}

		if (compressedObject == null) 
			return null;

		cachedObject = (T) ClassCompressor.byteToObject(compressedObject);
		return cachedObject;
	}

	public CompressedObject()
	{
	}

	public void beforeSerialize()
	{
		//		if (compressedObject != null && cachedObject != null)
		//			if (!OCHEMUtils.getMD5(compressedObject).equals(OCHEMUtils.getMD5(ClassCompressor.objectToByte((Serializable)cachedObject))))
		//				throw new RuntimeException("Gotcha");
		//				
		//		if (compressedObject == null && cachedObject != null)
		//			set(cachedObject);
	}

	public CompressedObject(T obj)
	{
		set(obj);
	}

	public static void main(String[] args)
	{
	}

}
