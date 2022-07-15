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

import java.io.IOException;
import java.sql.Timestamp;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.bson.types.ObjectId;

import qspr.workflow.utils.QSPRConstants;
import qspr.metaserver.transport.MongoConnection;

import com.mongodb.AggregationOptions;
import com.mongodb.BasicDBObject;
import com.mongodb.Cursor;
import com.mongodb.DBCollection;
import com.mongodb.DBCursor;
import com.mongodb.DBObject;
import com.mongodb.WriteConcern;

/**
 * An implementation of a MongoDB-based descriptors storage
 * @author novserj
 */
public class DescriptorsRepositoryImpl implements DescriptorsRepository
{
	private static Logger logger = LogManager.getLogger(DescriptorsRepositoryImpl.class);
	protected static String database = QSPRConstants.DESCRIPTORSCACHE;
	protected static String config_collection = "config";
	protected static String entry_collection = "entry";

	protected MongoConnection connection = null;

	final static int BATCH_TO_SAVE = 100;

	private String[] cutNulls(String[] identifiers)
	{
		if (identifiers == null)
			return null;
		int nonNullIdentifiers = 0;

		for (String identifier : identifiers)
			if (identifier != null)
				nonNullIdentifiers++;

		String[] newIdentifiers = new String[nonNullIdentifiers];
		nonNullIdentifiers = 0;

		for (String identifier : identifiers)
			if (identifier != null)
			{
				newIdentifiers[nonNullIdentifiers] = identifier;
				nonNullIdentifiers++;
			}
		return newIdentifiers;
	}

	private List<CacheEntry> getDescriptors(String[] identifiers, String fieldName, DescriptorConfigEntry config) throws Exception
	{
		DBCollection coll = connection.db().getCollection(entry_collection);

		BasicDBObject query = new BasicDBObject(); 
		query.append("config", new ObjectId(config.objectID));

		if (identifiers != null)
			query.append(fieldName, new BasicDBObject("$in", cutNulls(identifiers)));

		query.append("user", config.user);
		DBCursor cursor = coll.find(query);

		Map<String, CacheEntry> map = new HashMap<String, CacheEntry>();

		while (cursor.hasNext())
		{
			DBObject object = cursor.next();
			String identifier = object.get(fieldName).toString();
			CacheEntry entry = DBObjectHandler.getEntryFromDBObject(object, config);
			map.put(identifier, entry);
		}

		//Restore query order and fill in non-found items with nulls, if necessary
		List<CacheEntry> result = new ArrayList<CacheEntry>();
		if (identifiers != null)
		{
			for (String identifier : identifiers) 
			{
				if (identifier == null || map.get(identifier) == null)
					result.add(null);
				else
					result.add(map.get(identifier));
			}
		} else
		{
			List<String> ids = new ArrayList<String>();
			ids.addAll(map.keySet());
			Collections.sort(ids);
			for (String identifier : ids)
				result.add(map.get(identifier));
		}

		return result;		
	}

	public void updateDescriptors(DescriptorConfigEntry oldConfig,  DescriptorConfigEntry newConfig) throws Exception
	{
		DBCollection coll = connection.db().getCollection(entry_collection);

		saveConfig(oldConfig); // getting IDs in the database
		saveConfig(newConfig);

		BasicDBObject query = new BasicDBObject(); 
		query.append("config", new ObjectId(oldConfig.objectID));

		query.append("user", oldConfig.user);
		DBCursor cursor = coll.find(query);

		DBCollection coll1 = connection.db().getCollection(entry_collection);
		coll1.setWriteConcern(WriteConcern.UNACKNOWLEDGED);  // To speed up; even if some descriptors are not saved, this is not critical

		Set<String> cached = new HashSet<String>();

		int updated = 0, skipped = 0;

		while (cursor.hasNext())
		{
			DBObject object = cursor.next();
			CacheEntry entry = DBObjectHandler.getEntryFromDBObject(object, oldConfig); // object found, convert to entry

			if(cached.contains(entry.moleculeMD5)) {
				System.out.println("ignoring previously saved entry ");
				continue;
			}

			cached.add(entry.moleculeMD5);

			entry.config = newConfig; // update Cfg
			entry.objectID = null;
			entry.user = newConfig.user;

			DBObject obj = DBObjectHandler.getDBObjectFromEntry(entry);
			coll.insert(obj);
			updated++;
		}

		System.out.println("\nupdated="+updated + " skipped: " + skipped);
	}

	@Override
	synchronized public List<CacheEntry> getDescriptors(DescriptorConfigEntry config) throws Exception
	{
		return getDescriptors(null, "_id", config);
	}

	@Override
	synchronized public List<CacheEntry> getDescriptors(String[] md5, DescriptorConfigEntry config) throws Exception
	{
		return getDescriptors(md5, "md5", config);
	}

	@Override
	synchronized public List<CacheEntry> getDescriptors(Integer[] mp2, DescriptorConfigEntry config) throws Exception
	{
		String[] mp2string = new String[mp2.length];
		for (int i=0; i<mp2.length; i++)
			if (mp2[i] != null)
				mp2string[i] = mp2[i].toString();

		return getDescriptors(mp2string, "mp2", config);
	}

	@Override
	synchronized public void saveDescriptors(List<CacheEntry> cacheEntries) throws Exception
	{
		long timer = System.nanoTime();
		logger.info("Saving " + cacheEntries.size() + " new cache entries");
		DBCollection coll = connection.db().getCollection(entry_collection);

		coll.setWriteConcern(WriteConcern.UNACKNOWLEDGED);  // To speed up; even if some descriptors are not saved, this is not critical

		int i = 0, updated = 0;
		for (CacheEntry cacheEntry : cacheEntries) 
		{
			DBObject obj = DBObjectHandler.getDBObjectFromEntry(cacheEntry);
			if (cacheEntry.objectID == null)
				coll.insert(obj);
			else
			{
				obj.put("_id", new ObjectId(cacheEntry.objectID));
				coll.save(obj);// TODO: To Test
				updated++;
			}
			cacheEntry.objectID = ((ObjectId)obj.get("_id")).toString();

			if (++i % BATCH_TO_SAVE == 0)
				logger.info("Saved " + i + " new cache entries out of " + cacheEntries.size() + (updated>0? " including updated: "+updated:""));
		}
		timer = System.nanoTime() - timer;
		logger.info("Saved " + cacheEntries.size() + " entries in "+(timer/1000000)+"ms, an average of "+(timer / 1000000 / cacheEntries.size())+"ms per entry");
	}

	@Override
	synchronized public void saveConfig(DescriptorConfigEntry config) throws Exception
	{
		DescriptorConfigEntry entry = getDescriptorConfig(config.md5);

		if (entry != null)
		{
			config.objectID = entry.objectID;
			return;
		} else
		{
			DBCollection coll = connection.db().getCollection(config_collection);
			DBObject obj = DBObjectHandler.getDBObjectFromConfig(config);
			coll.insert(obj);
			config.objectID = obj.get("_id").toString();
		}
	}

	public void saveConfigWithMerge(DescriptorConfigEntry config) throws Exception
	{
		if (config.objectID == null) // Merge only makes sense for the config already in DB
		{
			saveConfig(config);
			return;
		}

		DescriptorConfigEntry entry = getDescriptorConfig(config.md5);

		if (entry == null) //No config with such md5... we just save/overwrite our current config
		{
			DBCollection coll = connection.db().getCollection(config_collection);
			DBObject obj = DBObjectHandler.getDBObjectFromConfig(config);
			obj.put("_id", new ObjectId(config.objectID));
			coll.save(obj);
			logger.info("Saved updated config "+config);
			return;
		} else
			if (!entry.objectID.equals(config.objectID)) // We have a config and it's not us... we need to merge entries to existing config and delete our config
			{

				DBCollection configColl = connection.db().getCollection(config_collection);
				DBCollection entryColl = connection.db().getCollection(entry_collection);

				DBObject query = new BasicDBObject("config", new ObjectId(config.objectID));
				DBObject update = new BasicDBObject("$set", new BasicDBObject("config", new ObjectId(entry.objectID)));
				entryColl.updateMulti(query, update);

				DBObject cquery = new BasicDBObject("id", new ObjectId(config.objectID));
				configColl.remove(cquery);

				logger.info("Merged "+config+" to "+entry+" due to md5 collision");
				config.objectID = entry.objectID;
			}
	}

	@Override
	synchronized public DescriptorConfigEntry getDescriptorConfig(String md5) throws Exception
	{
		DBCollection coll = connection.db().getCollection(config_collection) ;

		BasicDBObject query = new BasicDBObject(); 
		query.append("md5", md5);

		DBObject object = coll.findOne(query);
		return (object == null) ? null : DBObjectHandler.getConfigFromDBObject(object);
	}

	@Override
	synchronized public DescriptorConfigEntry getDescriptorConfigById(String objectId) throws Exception
	{
		DBCollection coll = connection.db().getCollection(config_collection) ;

		BasicDBObject query = new BasicDBObject(); 
		query.append("_id", new ObjectId(objectId));

		DBObject object = coll.findOne(query);
		return (object == null) ? null : DBObjectHandler.getConfigFromDBObject(object);
	}

	@Override
	public List<DescriptorConfigEntry> getAllConfigurations() throws Exception {
		return getConfigurations(null, true);
	}


	@Override
	public List<DescriptorConfigEntry> getConfigurations(String user) throws Exception
	{
		return getConfigurations(user, false);
	}

	synchronized private List<DescriptorConfigEntry> getConfigurations(String user, boolean all) throws Exception
	{
		long timer = System.nanoTime();
		DBCollection collConfig = connection.db().getCollection(config_collection) ;
		DBCollection collEntry = connection.db().getCollection(entry_collection) ;

		//clearLostEntries();

		List<DBObject> pipeline = new ArrayList<DBObject>();
		if(!all) pipeline.add(new BasicDBObject("$match", new BasicDBObject("user", user)));
		pipeline.add(new BasicDBObject("$project", new BasicDBObject("config", 1)
				.append("names_size", 1).append("values_size", 1)));
		pipeline.add(new BasicDBObject("$group", new BasicDBObject( "_id", "$config")
				.append("entries", new BasicDBObject("$sum", 1))
				.append("names_size", new BasicDBObject("$sum", "$names_size"))
				.append("values_size", new BasicDBObject("$sum", "$values_size"))));

		AggregationOptions options = AggregationOptions.builder().build();
		Cursor cursor = collEntry.aggregate(pipeline,options); 

		List<DescriptorConfigEntry> results = new ArrayList<DescriptorConfigEntry>();
		while (cursor.hasNext() )
		{
			DBObject res = cursor.next();
			DescriptorConfigEntry entry = DBObjectHandler.getConfigFromDBObject(collConfig.findOne(new BasicDBObject("_id", res.get("_id"))));
			entry.entriesCount = Long.valueOf(res.get("entries")+"");
			entry.entriesSize =  Long.valueOf(res.get("names_size")+"")+Long.valueOf(res.get("values_size")+"");
			results.add(entry);
		}

		logger.info("Got configurations for user "+user+" in "+(System.nanoTime() - timer) / 1000000 +"ms");
		return results;
	}

	@Override
	synchronized public void clearCache(String descType) throws Exception
	{
		List<DescriptorConfigEntry> configurations = getConfigurationsByType(descType);

		DBCollection collEntry = connection.db().getCollection(entry_collection);

		BasicDBObject query = null;
		for (DescriptorConfigEntry config : configurations) 
		{
			long timer = System.nanoTime();
			query = new BasicDBObject("config", new ObjectId(config.objectID));
			collEntry.remove(query);

			DBCollection collConfig = connection.db().getCollection(config_collection);
			query = new BasicDBObject("_id", new ObjectId(config.objectID));
			collConfig.remove(query);

			logger.info("Cleared cache entries by type "+descType+" for config "+config.objectID+" in "+(System.nanoTime() - timer) / 1000000 +"ms");
		}
	}

	@Override
	synchronized public void clearCache(String descType, String user) throws Exception
	{
		List<DescriptorConfigEntry> configurations = getConfigurationsByType(descType);

		DBCollection collEntry = connection.db().getCollection(entry_collection);

		BasicDBObject query = null;
		for (DescriptorConfigEntry config : configurations) 
		{
			long timer = System.nanoTime();
			query = new BasicDBObject("config", new ObjectId(config.objectID));
			query.append("user", user);
			collEntry.remove(query);

			query = new BasicDBObject("config", new ObjectId(config.objectID));
			long remainingEntries = collEntry.count(query);
			if (remainingEntries == 0)
			{
				DBCollection collConfig = connection.db().getCollection(config_collection);
				query = new BasicDBObject("_id", new ObjectId(config.objectID));
				collConfig.remove(query);
			}

			logger.info("Cleared cache entries by type "+descType+" for config "+config.objectID+" in "+(System.nanoTime() - timer) / 1000000 +"ms");
		}
	}

	@Override
	synchronized public void clearCache(DescriptorConfigEntry config) throws Exception
	{
		try{

			DBCollection collEntry = connection.db().getCollection(entry_collection);

			long timer = System.nanoTime();
			BasicDBObject query = new BasicDBObject("config", new ObjectId(config.objectID));
			//TODO remove once all caches are updated 
			if(config.user != null) query.append("user", config.user);  //If this is super user, it can delete all entries of the given type
			long initialEntries = collEntry.count(query);
			collEntry.remove(query);

			query = new BasicDBObject("config", new ObjectId(config.objectID));
			long remainingEntries = collEntry.count(query);
			if (remainingEntries == 0)
			{
				DBCollection collConfig = connection.db().getCollection(config_collection);
				query = new BasicDBObject("_id", new ObjectId(config.objectID));
				collConfig.remove(query);
			}

			logger.info("Cleared cached entries: " + initialEntries + " still left: " + remainingEntries + " by type " + config.type + 
					" for config "+config.objectID+" and user: " + (config.user == null ?"all" : config.user) + " in "+(System.nanoTime() - timer) / 1000000 +"ms");
		}catch(Exception ee){
			clearLostEntries();
			throw ee;
		}
	}

	@Override
	synchronized public List<DescriptorConfigEntry> getConfigurationsByType(String descType) throws Exception
	{
		long timer = System.nanoTime();
		List<DescriptorConfigEntry> entries = new ArrayList<DescriptorConfigEntry>();
		DBCollection coll = connection.db().getCollection(config_collection) ;

		BasicDBObject query = new BasicDBObject(); 
		query.append("type", descType);

		DBCursor cursor = coll.find(query);

		while (cursor.hasNext())
		{
			DBObject object = cursor.next();
			entries.add(DBObjectHandler.getConfigFromDBObject(object));
		}
		logger.info("Got "+entries.size()+" configurations by type "+descType+" in "+(System.nanoTime() - timer) / 1000000 +"ms");
		return entries;
	}

	private void cleanOrphanConfigurations() throws Exception
	{
		long timer = System.nanoTime();
		DBCollection collConfig = connection.db().getCollection(config_collection) ;
		DBCollection collEntry = connection.db().getCollection(entry_collection) ;

		List<DBObject> pipeline = new ArrayList<DBObject>();
		pipeline.add(new BasicDBObject("$project", new BasicDBObject("config", 1)));
		pipeline.add(new BasicDBObject("$group", new BasicDBObject( "_id", "$config")));

		//AggregationOutput output = collEntry.aggregate(pipeline);

		AggregationOptions options = AggregationOptions.builder().build();
		Cursor cursor = collEntry.aggregate(pipeline,options); 

		List<ObjectId> nonOrphans = new ArrayList<ObjectId>();
		while (cursor.hasNext() )
		{
			DBObject res = cursor.next();
			nonOrphans.add((ObjectId)res.get("_id"));
		}

		DBObject query =  new BasicDBObject("_id", new BasicDBObject("$nin", nonOrphans));
		collConfig.remove(query);
		logger.info("Cleaned orphan configurations in "+(System.nanoTime() - timer) / 1000000 +"ms");
	}

	private void cleanOrphanEntries() throws Exception
	{
		long timer = System.nanoTime();
		DBCollection collConfig = connection.db().getCollection(config_collection) ;
		DBCollection collEntry = connection.db().getCollection(entry_collection) ;

		List<DBObject> pipeline = new ArrayList<DBObject>();
		pipeline.add(new BasicDBObject("$project", new BasicDBObject("config", 1)));
		pipeline.add(new BasicDBObject("$group", new BasicDBObject( "_id", "$config")));

		AggregationOptions options = AggregationOptions.builder().build();
		Cursor cursor = collEntry.aggregate(pipeline,options); 

		logger.info("Aggregated orphan entries in "+(System.nanoTime() - timer) / 1000000 +"ms");

		while (cursor.hasNext() )
		{
			DBObject res = cursor.next();
			ObjectId id = (ObjectId)res.get("_id");
			DBObject config = collConfig.findOne(new BasicDBObject("_id", id));
			if (config == null)
				collEntry.remove(new BasicDBObject("config", id));
		}
		logger.info("Cleaned orphan entries in "+(System.nanoTime() - timer) / 1000000 +"ms");
	}

	@Override
	synchronized public void initCollections() throws Exception
	{
		DBCollection coll;

		coll = connection.db().getCollection(config_collection);
		coll.createIndex(new BasicDBObject("md5", 1), new BasicDBObject("unique", true));
		coll.createIndex(new BasicDBObject("type", 1));

		coll = connection.db().getCollection(entry_collection);
		coll.createIndex(new BasicDBObject("config", 1).append("md5", 1).append("user", 1));
		coll.createIndex(new BasicDBObject("config", 1).append("mp2", 1).append("user", 1));
		coll.createIndex(new BasicDBObject("user", 1));        
	}

	public DescriptorsRepositoryImpl()
	{
	}

	static class DBObjectHandler
	{
		private static CacheEntry getEntryFromDBObject(DBObject object, DescriptorConfigEntry config) throws IOException
		{
			CacheEntry entry = new CacheEntry();
			entry.objectID = ((ObjectId)object.get("_id")).toString();
			entry.attempts = (Integer)object.get("attempts");
			entry.error = (String)object.get("error");
			entry.config = config;
			entry.setMD5((String)object.get("md5"));

			if (object.get("mp2") != null)
				entry.mp2 = Integer.valueOf((String)object.get("mp2"));

			entry.user = (String)object.get("user");

			try{ //TODO fix parsing of time
				if (object.get("date_created") != null){
					String ss[] = object.get("date_created").toString().split("\\s+");
					DateFormat f = new SimpleDateFormat();
					entry.dateCreated = new Timestamp(f.parse(ss[0]).getTime());
					System.out.println(ss[0]);
				}
			}catch(ParseException e){

			}

			entry.zippedNames = (byte[])object.get("names");
			entry.zippedValues = (byte[])object.get("values");
			return entry;
		}

		private static DBObject getDBObjectFromEntry(CacheEntry entry) throws IOException
		{
			BasicDBObject object = new BasicDBObject();
			object.put("md5", entry.moleculeMD5);
			object.put("attempts", entry.attempts);
			object.put("error", entry.error);

			if (entry.mp2 != null)
				object.put("mp2", entry.mp2.toString());

			object.put("user", entry.user);
			object.put("date_created", entry.dateCreated.toString());
			object.put("config", new ObjectId(entry.config.objectID));


			if (entry.zippedNames == null || entry.zippedValues == null)
				entry.setNamesAndValues(new String[]{}, new float[]{}, false);

			object.put("names", entry.zippedNames);
			object.put("names_size", entry.zippedNames.length);

			object.put("values", entry.zippedValues);
			object.put("values_size", entry.zippedValues.length);
			return object;
		}

		private static DescriptorConfigEntry getConfigFromDBObject(DBObject object)
		{
			DescriptorConfigEntry config = new DescriptorConfigEntry();
			if(object == null){
				System.out.println("DBObject object is null!");
				return config;
			}
			config.objectID = (ObjectId)object.get("_id") == null? null : ((ObjectId)object.get("_id")).toString();
			config.md5 = (String)object.get("md5");
			config.description = (String)object.get("description");
			config.type = (String)object.get("type");
			config.setUser((String)object.get("user"));
			return config;
		}

		private static DBObject getDBObjectFromConfig(DescriptorConfigEntry config)
		{
			BasicDBObject object = new BasicDBObject();
			object.put("md5", config.md5);
			object.put("description", config.description);
			object.put("type", config.type);
			object.put("user", config.user);
			return object;
		}
	}

	@Override
	public void clearLostEntries() throws Exception {
		cleanOrphanConfigurations();
		cleanOrphanEntries();		
	}

}