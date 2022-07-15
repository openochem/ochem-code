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

package qspr.metaserver.transport;

import java.io.IOException;
import java.lang.reflect.InvocationHandler;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import qspr.util.OverloadControl;

public class ReattemptingProxy implements InvocationHandler
{
	Object obj;
	public ReattemptingProxy(Object obj)
	{ 
		this.obj = obj; 
	}
	
	public Object invoke(Object proxy, Method m, Object[] args) throws Throwable
	{
		try 
		{
			if (!m.isAnnotationPresent(Reattempt.class))
				return m.invoke(obj, args);
			Reattempt a = m.getAnnotation(Reattempt.class);
			OverloadControl oc = new OverloadControl(a.name(), a.timeout(), a.maxTimeout());
			while (true)
			{
				try
				{
			        return m.invoke(obj, args);
				} catch (InvocationTargetException e)
				{
					if (!(e.getTargetException() instanceof IOException))
						throw e;
					if (oc.getNumMaxTimeout() > a.numMaxTimeouts())
						throw e;
					oc.relax(e.getTargetException());
				}
			}
		} catch (InvocationTargetException e) 
		{
			throw e.getTargetException();
		} catch (Exception e) 
		{
			throw e;
		}
	}
}

