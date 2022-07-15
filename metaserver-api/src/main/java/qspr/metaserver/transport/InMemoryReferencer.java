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

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

/**
 * A simple in-memory data referencer
 * @author midnighter
 *
 */
public class InMemoryReferencer implements DataReferencer {

	/**
	 * The actual in-memory storage
	 */
	private static Map<String, Serializable> storage = new HashMap<String, Serializable>();
	
	@Override
	public DataReference saveReference(Serializable object, String collection) {
		InMemoryReference reference = new InMemoryReference(UUID.randomUUID().toString());
		storage.put(reference.getReference(), object);
		
		return reference;
	}

	@Override
	public Serializable getReference(DataReference reference) throws DataReferenceException {
		InMemoryReference ref = (InMemoryReference) reference;
		if (!storage.containsKey(ref.getReference()))
			throw new DataReferenceException("Unknown in memory reference: " + reference);
		
		return storage.get(ref.getReference());
	}

	@Override
	public Long getDataSize(DataReference reference) {
		// TODO Auto-generated method stub
		return 0l;
	}

	@Override
	public DataReference putDataBytes(byte[] data, 
			String collection) throws DataReferenceException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public byte[] getDataBytes(DataReference reference)
			throws DataReferenceException {
		// TODO Auto-generated method stub
		return null;
	}

}
