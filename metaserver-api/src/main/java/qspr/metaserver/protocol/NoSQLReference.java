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

package qspr.metaserver.protocol;

import qspr.metaserver.transport.DataReference;
import qspr.metaserver.transport.NoSqlTransport;

public class NoSQLReference implements DataReference
{
	private static final long serialVersionUID = 1L;
	
	public String collection;
	public String database;
	private String key;
	
	public NoSQLReference(String key, String database)
	{
		this.key = key;
		this.database = database;
		collection = NoSqlTransport.DEFAULT_COLLECTION;
	}
	
	public NoSQLReference(String key, String database, String collection)
	{
		this.key = key;
		this.database = database;
		this.collection = collection;
	}
	
	public String toString()
	{
		return database+":"+collection+":"+key;
	}

	@Override
	public String getReference() {
		return key;
	}
}
