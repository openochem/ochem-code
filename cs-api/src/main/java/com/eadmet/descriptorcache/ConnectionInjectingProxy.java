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

import java.lang.reflect.InvocationHandler;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;

import qspr.metaserver.transport.MongoConnection;
import qspr.metaserver.transport.MongoConnectionFactory;
import qspr.metaserver.transport.NoSqlTransport;

public class ConnectionInjectingProxy implements InvocationHandler
{
	DescriptorsRepositoryImpl obj;

	public ConnectionInjectingProxy(DescriptorsRepositoryImpl obj)
	{ 
		this.obj = obj; 
	}

	public Object invoke(Object proxy, Method m, Object[] args) throws Throwable
	{
		MongoConnection c = null;
		try 
		{
			c = MongoConnectionFactory.getDbPersistentConnection(NoSqlTransport.host, DescriptorsRepositoryImpl.database);
			obj.connection = c;
			return m.invoke(obj, args);
		} catch (InvocationTargetException e) 
		{
			throw e.getTargetException();
		} catch (Exception e) 
		{
			throw e;
		} finally
		{
			if (c != null)
			{
				c.mongo().close();
				c.clear();
			}
		}
	}

	//	private static final Logger logger = Logger.getLogger(ConnectionInjectingProxy.class);
}

