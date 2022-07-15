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
import java.lang.reflect.Proxy;

public class MongoConnectionFactory
{
	/**
	 * Returns persistent connection
	 * @param host
	 * @param database
	 * @return
	 * @throws IOException
	 */
	public static MongoConnection getDbPersistentConnection(String host, String database) throws IOException
	{
		MongoConnection mpc = new MongoPersistentConnectionImpl();
		MongoConnection proxy = (MongoConnection)Proxy.newProxyInstance(mpc.getClass().getClassLoader(), new Class[]{MongoConnection.class}, new ReattemptingProxy(mpc));
		proxy.getConnection(host, database);
		return proxy;
	}

	/**
	 * Returns one time connection
	 * @param host
	 * @param database
	 * @return
	 * @throws IOException
	 */
	public static MongoConnection getDbNoProxyConnection(String host, String database) throws IOException
	{
		MongoConnection mpc = new MongoPersistentConnectionImpl();
		mpc.getOneTimeConnection(host, database);
		return mpc;
	}
}