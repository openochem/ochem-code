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

package com.eadmet.mmpa.domain;

/**
 * A "thin" representation of MoleculerMatchedPair class.
 * Use to efficiently keep pairs in memory.
 * 
 * @author midnighter
 *
 */
public class ThinPair
{
	public long id;
	public int mol1Id;
	public int mol2Id;
	public long transformationId;
	
	public void invert()
	{
		int tmp = mol1Id;
		mol1Id = mol2Id;
		mol2Id = tmp;
	}
	
	public ThinPair() {
		
	}

	public ThinPair(long id, int mol1Id, int mol2Id, long transformationId)
	{
		this.id = id;
		this.mol1Id = mol1Id;
		this.mol2Id = mol2Id;
		this.transformationId = transformationId;
	}
	
	public MMPair getPair() {
		return MMPair.thin(id, mol1Id, mol2Id);
	}
}
