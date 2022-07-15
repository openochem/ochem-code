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

package qspr.entities;


import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlValue;
import javax.xml.bind.annotation.adapters.XmlAdapter;
import javax.xml.bind.annotation.adapters.XmlJavaTypeAdapter;

@XmlRootElement
public class AnalogProfile {
	@XmlAttribute
	public String id = null;
	
	@XmlJavaTypeAdapter(TreeMapAdapter.class)
	public TreeMap<Long, String> profileTree = null;
	
	@XmlAttribute
	public float pearson = 0.0f;

	@XmlAttribute
	public float similarityCount = 0;
	
	public ModelMapping modelMapping;
	
	@XmlElement
	public List<String> profile = new ArrayList<String>();
	
	@Override
	public String toString() {
		return id + " - pearson: " + pearson;
	}
    
    // treemap marshalling
    public static class MapType {
        @XmlElement(name ="entry")
        public List<MapEntryTree> entryList = new ArrayList<MapEntryTree>();
    }
    
    public static class MapEntryTree {
        @XmlAttribute
        public Long key;
        @XmlValue
        public String value;
    }
    
    public static class TreeMapAdapter extends XmlAdapter<MapType, Map<Long, String>> {
        @Override
        public MapType marshal(Map<Long, String> map) {
            MapType mapType = new MapType();
            for (Entry<Long, String> entry : map.entrySet()) {
                MapEntryTree mapEntryTree = new MapEntryTree();
                mapEntryTree.key = entry.getKey();
                mapEntryTree.value = entry.getValue();
                mapType.entryList.add(mapEntryTree);
            }
            return mapType;
        }
        
        @Override
        public Map<Long, String> unmarshal(MapType type) throws Exception {
            Map<Long, String> map = new TreeMap<Long, String>();
            for (MapEntryTree entry : type.entryList) { 
                map.put(entry.key, entry.value);
            }
            return map;
        }
    }
    
//	public static class HashMapGenericAdapter<K,V> extends XmlAdapter<MapType<K,V>, Map<K, V>> {
//        @Override
//        public MapType<K,V> marshal(Map<K, V> map) {
//        	MapType<K,V> mapType = new MapType<K,V>();
//            for (Entry<K, V> entry : map.entrySet()) {
//            	MapEntry<K,V> mapEntry = new MapEntry<K,V>();
//                mapEntry.key = entry.getKey();
//                mapEntry.value = entry.getValue();
//                mapType.entryList.add(mapEntry);
//            }
//            return mapType;
//        }
//        @Override
//        public Map<K, V> unmarshal(MapType<K,V> type) throws Exception {
//            Map<K, V> map = new HashMap<K, V>();
//            for (MapEntry<K,V> entry : type.entryList) { 
//                map.put(entry.key, entry.value);
//            }
//            return map;
//        }
//    }
// 
//    public static class MapType<K,V> {
//        @XmlElement(name ="entry")
//        public List<MapEntry<K,V>> entryList = new ArrayList<MapEntry<K,V>>();
//    }
// 
//    public static class MapEntry<K,V> {
//        @XmlAttribute
//        public K key;
//        @XmlValue
//        public V value;
//    }
//    
}

