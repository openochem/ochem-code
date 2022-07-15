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

package qspr.util;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringEscapeUtils;

public class StringUtils 
{
	public static String listToString(List<String> list)
	{
		StringBuffer sb = new StringBuffer();
		for (String s : list)
			sb.append(StringEscapeUtils.escapeCsv(StringEscapeUtils.escapeJava(s))+",");
		return org.apache.commons.lang3.StringUtils.chomp(sb.toString(), ",");
	}
	
	public static String arrayToString(String[] list)
	{
		StringBuffer sb = new StringBuffer();
		for (String s : list)
			sb.append(StringEscapeUtils.escapeCsv(StringEscapeUtils.escapeJava(s))+",");
		return org.apache.commons.lang3.StringUtils.chomp(sb.toString(), ",");
	}
	
	public static List<String> stringToList(String data)
	{
		List<String> result = new ArrayList<String>();
		if (data == null || data.equals(""))
			return result;
		String[] pieces = data.split(",(?=([^\"]*\"[^\"]*\")*[^\"]*$)");
		for (String piece : pieces)
		{
			String value = StringEscapeUtils.unescapeJava(StringEscapeUtils.unescapeCsv(piece));
			if (value.equals("null"))
				result.add(null);
			else
				result.add(value);
		}
		if (data.endsWith(","))
			result.add("");
		return result;
	}
	
	public static String mapToString(Map<String, String> map)
	{
		StringBuffer sb = new StringBuffer();
		for (String s : map.keySet())
			sb.append(s+":"+StringEscapeUtils.escapeCsv(StringEscapeUtils.escapeJava(map.get(s)))+",");
		return org.apache.commons.lang3.StringUtils.chomp(sb.toString(), ",");
	}
	
	public static Map<String,String> stringToMap(String data)
	{
		Map<String,String> result = new HashMap<String,String>();
		String[] pieces = data.split(",(?=([^\"]*\"[^\"]*\")*[^\"]*$)");
		for (String piece : pieces)
		{
			String[] keyvalue = piece.split(":",2);
			if (keyvalue[1].equals("null"))
				keyvalue[1] = null;
			result.put(keyvalue[0], StringEscapeUtils.unescapeJava(StringEscapeUtils.unescapeCsv(keyvalue[1])));
		}
		return result;
	}
}
