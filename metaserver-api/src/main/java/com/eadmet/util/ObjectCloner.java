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

package com.eadmet.util;

import java.io.*;
public class ObjectCloner
{
	// so that nobody can accidentally create an ObjectCloner object
	private ObjectCloner(){}
	// returns a deep copy of an object
	static public Object deepCopy(Object oldObj) throws Exception
	{
		ObjectOutputStream oos = null;
		ObjectInputStream ois = null;
		try
		{
			ByteArrayOutputStream bos = 
					new ByteArrayOutputStream(); // A
					oos = new ObjectOutputStream(bos); // B
					// serialize and pass the object
					oos.writeObject(oldObj);   // C
					oos.flush();               // D
					ByteArrayInputStream bin = 
							new ByteArrayInputStream(bos.toByteArray()); // E
					ois = new ObjectInputStream(bin);                  // F
					// return the new object
					return ois.readObject(); // G
		}
		catch(Exception e)
		{
			System.out.println("Exception in ObjectCloner = " + e);
			throw(e);
		}
		finally
		{
			oos.close();
			ois.close();
		}
	}

}