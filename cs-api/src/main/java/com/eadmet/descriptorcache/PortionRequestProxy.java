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

import java.lang.reflect.Array;
import java.lang.reflect.InvocationHandler;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.List;

public class PortionRequestProxy implements InvocationHandler
{
	Object obj;
	public PortionRequestProxy(Object obj)
	{ 
		this.obj = obj; 
	}
	
	public Object invoke(Object proxy, Method m, Object[] args) throws Throwable
	{
		try 
		{
			if (!m.isAnnotationPresent(Portion.class))
				return m.invoke(obj, args);
			Portion a = m.getAnnotation(Portion.class);
			Object o = args[0];
			Class<?> c = o.getClass();
	
			if (!c.isArray())
				throw new Exception("Portion requests support only array arguments (can be extended to lists if necessary)");
			
			if (!List.class.isAssignableFrom(m.getReturnType()))
				throw new Exception("Portion requests support only list return types");
			
			int length = Array.getLength(o);
			List<Object> result = new ArrayList<Object>();
			int blockSize = a.size();
			int blockCount = length / blockSize + ((length % blockSize == 0) ? 0 : 1);
			for (int block = 0; block < blockCount; block++)
			{
				int curSize = Math.min(length - blockSize * block, blockSize);
//				System.out.println("Requesting portion "+block*blockSize+" - "+(block*blockSize + curSize));
				Object subarg = Array.newInstance(c.getComponentType(), curSize);
				for (int i = 0; i < curSize; i++)
					Array.set(subarg, i, Array.get(o, block * blockSize + i));
				
				args[0] = subarg;
				List<?> subres = (List<?>)m.invoke(obj, args);
				result.addAll(subres);
			}
			args[0] = o;
			return result;
		} catch (InvocationTargetException e) 
		{
			throw e.getTargetException();
		} catch (Exception e) 
		{
			throw e;
		}
	}
}

