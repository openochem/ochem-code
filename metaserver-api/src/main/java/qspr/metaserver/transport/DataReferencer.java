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

/**
 * An abstract data referencer/dereferencer.
 * Can:
 * - store Serializable objects in abstract storage, 
 * - give back an abstract reference
 * - dereference objects using the data reference
 * 
 * @author midnighter
 */
public interface DataReferencer
{
	public DataReference putDataBytes(byte [] data, String collection) throws DataReferenceException;
	public byte[] getDataBytes(DataReference reference) throws DataReferenceException;
	public DataReference saveReference(Serializable object, String collection) throws DataReferenceException;
	public Serializable getReference(DataReference reference) throws DataReferenceException;
	
	/**
	 * Return size in bytes or null if reference does not exist
	 * @param reference
	 * @return
	 */
	
	public Long getDataSize(DataReference reference);
}
