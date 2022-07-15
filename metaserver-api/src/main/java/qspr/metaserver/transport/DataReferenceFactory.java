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

import com.eadmet.exceptions.UserFriendlyException;

import qspr.metaserver.protocol.NoSQLReference;

public class DataReferenceFactory {
	public static Class<? extends DataReferencer> defaultReferencer = NoSqlTransport.class;

	public static DataReferencer createReferencer(DataReferencerType type) {
		try{
			switch(type){
			case MEMORY:
				return InMemoryReferencer.class.newInstance();
			case NOSQL:
				return NoSqlTransport.class.newInstance();
			}
		}catch (Exception e){
			throw new UserFriendlyException(e.getMessage());
		}
		throw new UserFriendlyException("type is not initialised" + type);
	}
	
	public static DataReferencer createReferencer() {
		try{
			return defaultReferencer.newInstance();
		}catch (Exception e){
			throw new UserFriendlyException(e.getMessage());
		}
	}

	/**
	 * Retrieves reference from the database
	 * @param key
	 * @param database
	 * @return
	 * @throws DataReferenceException in case if reference does not exist or there is an error
	 */
	
	public static DataReference resolveReference(String key, String database) throws DataReferenceException {
		DataReference ref= createReference(key, database);
		DataReferencer referencer = createReferencer();
		Long size = referencer.getDataSize(ref);
		if(size == null)throw new DataReferenceException("Reference does not exist");
		return ref;
	}

	/**
	 * Creates a reference that can be used to store or check data
	 * @param key
	 * @param database
	 * @return
	 * @throws DataReferenceException
	 */
	
	public static DataReference createReference(String key, String database) throws DataReferenceException {
			DataReferencer ref = createReferencer();
			if(ref instanceof NoSqlTransport)
				return new NoSQLReference(key, database);
			else
				if(ref instanceof InMemoryReferencer)
					return new InMemoryReference(key);
				else
					throw new DataReferenceException("Unknown referencer type: " + ref);
	}
	
	public static DataReferencer createReferencer(DataReference ref) throws DataReferenceException {
		if (ref instanceof InMemoryReference)
			return new InMemoryReferencer();
		else if (ref instanceof NoSQLReference)
			return new NoSqlTransport();
		else
			throw new DataReferenceException("Unknown reference type: " + ref);
	}
	
	public enum DataReferencerType {
		NOSQL, MEMORY
	}
	
}
