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

import com.mongodb.DB;
import com.mongodb.MongoClient;
import com.mongodb.MongoClientOptions;
import com.mongodb.MongoClientURI;
import com.mongodb.ServerAddress;
import com.mongodb.WriteConcern;

class MongoPersistentConnectionImpl implements MongoConnection
{
	MongoClient mongo;
	DB db;
	//MongoDatabase basa;

	@Override
	public void getConnection(String host, String database) throws IOException
	{
		MongoClientOptions options = 
				MongoClientOptions.builder()
				.connectionsPerHost(3)
				.connectTimeout(120000)
				.build();

		MongoClientURI uri = new MongoClientURI(host+"/?authSource="+database);
		mongo = new MongoClient(new ServerAddress(uri.getHosts().get(0)),uri.getCredentials(), options);
		db = mongo.getDB(database);
		db.setWriteConcern(WriteConcern.JOURNALED);
		//basa = mongo.getDatabase(database).withWriteConcern(WriteConcern.JOURNALED);
	}

	@Override
	public void getOneTimeConnection(String host, String database) throws IOException {
		MongoClientURI uri = new MongoClientURI(host+"/?authSource="+database);
		mongo = new MongoClient(uri);
		db = mongo.getDB(database);
		//basa = mongo.getDatabase(database);
	}

	@Override
	public MongoClient mongo() 
	{
		return mongo;
	}

	@Override
	public DB db() 
	{
		return db;
	}

	@Override
	public void clear() 
	{
		mongo = null;
		db = null;
	}

	public static void main(String[] args) throws Exception{
		String hosts [] = {"mongodb://tom:jerry@eco", "mongodb://cpu"};

		String database ="descriptorscache";

		for (String host:hosts) {

			MongoConnection c  = MongoConnectionFactory.getDbNoProxyConnection(host, database);
			for (String collection : c.db().getCollectionNames()) 
			{
				System.out.println(host+ ": "+ collection);
			}
			c.mongo().close();

			c  = MongoConnectionFactory.getDbPersistentConnection(host, database);
			for (String collection : c.db().getCollectionNames()) 
			{
				System.out.println(host+ ": "+ collection);
			}
		}
	}

}