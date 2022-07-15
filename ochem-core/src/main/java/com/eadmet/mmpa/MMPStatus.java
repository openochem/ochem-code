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

package com.eadmet.mmpa;

import java.util.LinkedHashMap;
import java.util.Map;

import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement
public class MMPStatus
{
	public Map<String, Long> molCount = new LinkedHashMap<String, Long>();
	public Map<String, TableInfo> tableInfo = new LinkedHashMap<String, TableInfo>();
	
	@XmlElement
	public long totalSize = 0;
	
	@XmlElement
	public int getUnindexedMolecules() {
		int count = 0;
		for (String key : molCount.keySet())
			if (!"Indexed".equals(key))
				count += molCount.get(key);
		return count;
	}
	
	public TableInfo getTableInfo(String tableName) {
		TableInfo info = tableInfo.get(tableName);
		if (info == null)
			tableInfo.put(tableName, info = new TableInfo());
		
		return info;
	}
	
	public boolean indexingJobRunning;
}

class TableInfo {
	@XmlElement
	long size;
	
	@XmlElement
	long dataSize;
	
	@XmlElement
	long indexSize;
	
	@XmlElement
	long count;
}
