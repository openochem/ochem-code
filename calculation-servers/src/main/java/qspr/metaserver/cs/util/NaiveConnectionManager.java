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

package qspr.metaserver.cs.util;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;

public class NaiveConnectionManager
{
	Connection conn;
	public String dbUrl = null;
	public String dbPassword = null;
	public String dbUsername = null;

	Map<String, PreparedStatement> statements = null;

	private boolean isConnectionActive() throws SQLException
	{
		return (conn != null && !conn.isClosed());
	}

	private Connection getConnection() throws SQLException, InterruptedException
	{
		if (isConnectionActive())
			return conn;
		return conn = DriverManager.getConnection(dbUrl, dbUsername, dbPassword);
	}

	public void closeConnection() throws SQLException
	{
		if (isConnectionActive())
			conn.close();
	}

	public PreparedStatement getStatement(String sql, boolean cache) throws SQLException, InterruptedException {
		if(!isConnectionActive()) {
			statements = new HashMap<String,PreparedStatement>();
			getConnection();
		}

		if(!isConnectionActive())return null;

		if(!cache) return conn.prepareStatement(sql);;

		if(!statements.containsKey(sql))
			statements.put(sql,conn.prepareStatement(sql));
		return statements.get(sql);
	}
}

