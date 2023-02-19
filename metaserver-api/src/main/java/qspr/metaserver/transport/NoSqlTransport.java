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

import java.io.BufferedInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.bson.types.ObjectId;

import qspr.metaserver.protocol.NoSQLReference;
import qspr.util.ClassCompressor;

import com.eadmet.utils.OCHEMUtils;
import com.eadmet.utils.config.ConfigurableClass;
import com.eadmet.utils.config.ConfigurableProperty;
import com.mongodb.BasicDBObject;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import com.mongodb.gridfs.GridFS;
import com.mongodb.gridfs.GridFSDBFile;
import com.mongodb.gridfs.GridFSInputFile;

/**
 * A MongoDB implementation of the key-value storage transport.
 * 
 *  TODO Sergii: Think of what functionality should be exposed via an abstract transport interface.
 *  Currently there is a lot of MongoDB-specific functionality / Midnighter
 *  
 * @author novserj
 *
 */
@ConfigurableClass(name = "mongodb", comment = "MongoDB warehouse connectivity options")
public class NoSqlTransport implements DataReferencer
{
	private static transient Logger logger = LogManager.getLogger(NoSqlTransport.class);

	public static final String DESCRIPTORSCACHE = "descriptorscache";

	@ConfigurableProperty(name = "host")
	public static String host = "undefined";

	@ConfigurableProperty(name = "username")
	public static String username = "root";

	@ConfigurableProperty(name = "password")
	public static String password = "";

	public static int trials = 12;

	public final static String DEFAULT_COLLECTION ="Data";

	public final static String MONGOID ="_id";

	static 
	{
		System.setProperty("socksNonProxyHosts", "*");
	}


	/** For future changes, some useful commands
	 * 		MongoDatabase db = mongo.getDatabase(reference.database);
		GridFSBucket fs = GridFSBuckets.create(db, reference.collection);

		final Document doc = new Document("myKey", "myValue");
		final String jsonString = doc.toJson();
		final Document doc = Document.parse(jsonString);

	 */
	private static void ensurePowerOf2Collections(String database) throws IOException
	{
		if (Math.random() > 0.1) // Execute it every one out of ten times... it's very fast, but non the less...
			return;

		MongoConnection c  = MongoConnectionFactory.getDbNoProxyConnection(host, database);

		System.out.println();
		for (String collection : c.db().getCollectionNames()) 
		{
			DBObject cmd = new BasicDBObject();
			cmd.put("collMod", collection);
			cmd.put("usePowerOf2Sizes", true);
			c.db().command(cmd);
		}
		c.mongo().close();
		c.clear();
	}

	private static Long getDataSize(NoSQLReference reference) throws IOException
	{
		MongoConnection c  = MongoConnectionFactory.getDbNoProxyConnection(host, reference.database);
		GridFS fs = new GridFS(c.db(), reference.collection);
		GridFSDBFile out = fs.findOne(new BasicDBObject( MONGOID , new ObjectId(reference.getReference())));
		Long res = out == null? 0:out.getLength();
		c.mongo().close();
		c.clear();
		return res;
	}


	public static byte[] getDataSafely(NoSQLReference reference) throws IOException
	{
		int attempts = 0;
		while (true)
			try
		{
				return getData(reference);
		} catch (Exception e)
		{
			try {
				if (attempts++ < 12)
					Thread.sleep(1000*attempts);
				else
					throw new IOException(e);
			} catch (InterruptedException ie)
			{
				logger.error("Cannot resolve a MongoDB reference. Host is "+ host);
				throw new IOException(e);
			}
		}	
	}

	public byte[] getDataBytes(DataReference reference) throws DataReferenceException
	{
		try{
			return getDataSafely((NoSQLReference) reference);
		}catch(IOException e){
			throw new DataReferenceException(e.getMessage());
		} 
	}

	public NoSQLReference putDataBytes(byte[] data, String database)
	{
		try{
			return putDataSafely(data, database, DEFAULT_COLLECTION, null);
		}catch(IOException e){
			throw new DataReferenceException(e.getMessage());
		}
	}

	public static NoSQLReference putDataSafely(byte[] data, String database, String collection) throws IOException
	{
		return putDataSafely(data, database, collection, null);
	}

	public static NoSQLReference putDataSafely(byte[] data, String database, String collection, Map<String, Object> metadata) throws IOException
	{
		int attempts = 0;
		while (true)
			try
		{
				String md5 = OCHEMUtils.getMD5(data);
				NoSQLReference ref = findByMD5(database, collection, md5);
				if(ref != null) return ref;				
				return putData(data, database, collection, metadata, md5);
		} catch (Exception e)
		{
			try {
				if (attempts++ < trials)
					Thread.sleep(1000*attempts);
				else
					throw new IOException(e);
			} catch (InterruptedException ie)
			{
				throw new IOException(e);
			}
		}	
	}

	private static NoSQLReference putData(byte[] data, String database, String collection, Map<String, Object> metadata, String md5) throws Exception
	{
		long time = Calendar.getInstance().getTimeInMillis();

		MongoConnection c  = MongoConnectionFactory.getDbNoProxyConnection(host, database);
		GridFS fs = new GridFS(c.db(), collection);
		GridFSInputFile in = fs.createFile(data);
		in.put("md5", md5);

		if (metadata != null)
			for (String key : metadata.keySet())
				in.put(key, metadata.get(key));

		in.save();
		c.mongo().close();
		NoSQLReference reference = new NoSQLReference(in.getId().toString(), database, collection);
		ensurePowerOf2Collections(reference.database);

		logger.info("Saved "+data.length+" bytes to MongoDB["+reference+", "+md5+"] in " + (Calendar.getInstance().getTimeInMillis() - time) + "ms");
		return reference;
	}

	private static Map<String, String> mongoDbFileToMap(GridFSDBFile file)
	{
		if (file == null)
			return null;
		Map<String, String> res = new HashMap<String, String>();
		res.put("filename", file.getFilename());
		res.put("size", ""+file.getLength());
		res.put(MONGOID, file.getId().toString());

		if (file.get("task_id") != null)
			res.put("task_id", file.get("task_id").toString());

		if (file.get("parent_task_id") != null)
			res.put("parent_task_id", file.get("parent_task_id").toString());

		if (file.get("server_address") != null)
			res.put("server_address", file.get("server_address").toString());    	

		if (file.get("server_name") != null)
			res.put("server_name", file.get("server_name").toString());    	

		if (file.get("server_configuration") != null)
			res.put("server_configuration", file.get("server_configuration").toString());    	

		res.put("date", file.getUploadDate().toString());

		return res;
	}

	public static List<Map<String, String>> listFiles(String database, String collection) throws IOException
	{
		MongoConnection c  = MongoConnectionFactory.getDbNoProxyConnection(host, database);
		GridFS fs = new GridFS(c.db(), collection);
		List<GridFSDBFile> out = fs.find(new BasicDBObject());
		List<Map<String, String>> result = new ArrayList<Map<String, String>>();
		for (GridFSDBFile file : out)
			result.add(mongoDbFileToMap(file));
		c.mongo().close();
		return result;
	}

	public static Map<String, String> listFile(String database, String collection, String id) throws IOException
	{
		MongoConnection c  = MongoConnectionFactory.getDbNoProxyConnection(host, database);
		GridFS fs = new GridFS(c.db(), collection);
		GridFSDBFile out = fs.findOne(new BasicDBObject(MONGOID, new ObjectId(id)));
		c.mongo().close();
		return mongoDbFileToMap(out);
	}

	public static void deleteFile(String database, String collection, String id) throws IOException
	{
		MongoConnection c  = MongoConnectionFactory.getDbNoProxyConnection(host, database);
		GridFS fs = new GridFS(c.db(), collection);
		fs.remove(new ObjectId(id));
		c.mongo().close();
	}


	private static NoSQLReference findByMD5(String database, String collection, String md5) throws Exception
	{
		MongoConnection c  = MongoConnectionFactory.getDbNoProxyConnection(host, database);
		GridFS fs = new GridFS(c.db(), collection);
		GridFSDBFile out = fs.findOne(new BasicDBObject( "md5" ,md5));
		NoSQLReference reference = out == null? null: 
			new NoSQLReference(out.getId().toString(), database, collection);
		c.mongo().close();
		return reference;
	}

	private static byte[] getData(NoSQLReference reference) throws Exception
	{
		ensurePowerOf2Collections(reference.database);
		long time = Calendar.getInstance().getTimeInMillis();
		MongoConnection c  = MongoConnectionFactory.getDbNoProxyConnection(host, reference.database);
		GridFS fs = new GridFS(c.db(), reference.collection);
		GridFSDBFile out = null;

		out = fs.findOne(new BasicDBObject( MONGOID , new ObjectId(reference.getReference())));

		if (out == null){
			c.mongo().close();
			throw new IOException("Cannot resolve MongoDB reference " + reference);
		}

		if(out.getLength()>Integer.MAX_VALUE)
			throw new IOException("The length of MongoDB reference " + out.getLength() + " > Integer.MAX_VALUE; it can't be resolved");

		ByteArrayOutputStream baos = new ByteArrayOutputStream((int)out.getLength());
		BufferedInputStream bis = new BufferedInputStream(out.getInputStream());
		byte[] buffer = new byte[16536];
		int count = -1;
		while ((count = bis.read(buffer)) != -1)
			baos.write(buffer, 0, count);
		bis.close();
		baos.close();
		c.mongo().close();
		byte[] res = baos.toByteArray();
		logger.info("Downloaded "+res.length+" bytes from MongoDB ("+reference.collection+") in " + (Calendar.getInstance().getTimeInMillis() - time) + "ms.");
		return res;
	}


	public static Set<String> getReferenceIds(String database, String collection) throws Exception
	{
		MongoConnection c  = MongoConnectionFactory.getDbNoProxyConnection(host, database);
		DBCollection coll = c.db().getCollection(collection+".files") ;
		DBCursor cursor = coll.find(); // all references in the collection
		Set<String> ids = new HashSet<String>();

		while (cursor.hasNext())
		{
			DBObject object = cursor.next();
			ids.add(object.get(MONGOID).toString());
		}

		cursor.close();
		c.mongo().close();
		return ids;
	}

	public static void deleteReferenceByIds(String database, String collection, List<String> referenceIds)  throws IOException
	{
		long time = Calendar.getInstance().getTimeInMillis();
		MongoConnection c  = MongoConnectionFactory.getDbNoProxyConnection(host, database);
		GridFS fs = new GridFS(c.db(), collection);
		for(String id:referenceIds)
			fs.remove(new ObjectId(id));
		c.mongo().close();
		logger.info("Deleted "+referenceIds.size()+" records in " + (Calendar.getInstance().getTimeInMillis() - time) + "ms.");
	}

	@Override
	public DataReference saveReference(Serializable object, String database) throws DataReferenceException {
		byte[] bytes = ClassCompressor.objectToByte(object);
		try {
			return putDataSafely(bytes, database, DEFAULT_COLLECTION, null);
		} catch (IOException e) {
			throw new DataReferenceException(e);
		}
	}

	@Override
	public Serializable getReference(DataReference reference)
			throws DataReferenceException {
		try {
			NoSQLReference ref = (NoSQLReference) reference;
			if(ref.collection == null) ref.collection = DEFAULT_COLLECTION;
			byte[] bytes = getDataSafely(ref);
			return ClassCompressor.byteToObject(bytes);
		} catch (IOException e) {
			throw new DataReferenceException(e);
		}
	}

	@Override
	public Long getDataSize(DataReference reference){
		try {
			return getDataSize((NoSQLReference)reference);
		} catch (IOException e) {
			return null;
		}
	}

	static public void specificCacheCleaningjobTest() throws Exception
	{

		MongoConnection c  = MongoConnectionFactory.getDbPersistentConnection("qspr.helmholtz-muenchen.de", DESCRIPTORSCACHE);
		DBCollection collConfig = c.db().getCollection("config") ;
		DBObject configQuery = new BasicDBObject("type", "CDK");
		List<ObjectId> configs = new ArrayList<ObjectId>();
		DBCursor res = collConfig.find(configQuery);
		for (DBObject result : res) 
		{
			System.out.println("Found ObjectID of interest: "+(ObjectId)result.get(MONGOID));
			configs.add((ObjectId)result.get(MONGOID));
		}

		DBCollection collEntry = c.db().getCollection("entry");

		for (ObjectId oid : configs) 
		{
			System.out.println("Processing "+oid);
			DBObject entryQuery = new BasicDBObject("config", oid).append("values_size", 24).append("error", null);
			DBObject fields = new BasicDBObject("md5", 1).append("mp2", 1).append("names_size", 1).append("values_size", 1);
			DBCursor eres = collEntry.find(entryQuery, fields);
			for (DBObject result : eres) 
			{
				System.out.println(result.get("md5")+" "+result.get("mp2")+" "+result.get("names_size")+" "+result.get("values_size"));
				collEntry.remove(result);
			}
		}
		c.mongo().close();
	}


	public static void main(String[] args) throws IOException 
	{
		NoSqlTransport.host = "qsar";

		byte data[] = "552bc3dbe4b05787bacd3cb6".getBytes();

		System.out.println(OCHEMUtils.getMD5(data));

		NoSQLReference ref = putDataSafely(data, "metaserver","Data");
		System.out.println(ref);

		ref = putDataSafely(data, "metaserver","Data");
		System.out.println(ref);

		List<String> s = new ArrayList<String>();
		s.add(ref.getReference());
		deleteReferenceByIds("metaserver","Data", s);

		ref = putDataSafely(data, "metaserver","Data");
		System.out.println(ref);
	}

}
