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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import qspr.entities.Property;

public class SmartText 
{
	public static String transform(String text)
	{
		Pattern p = Pattern.compile("PROP\\$([0-9]+)");
		Matcher m = p.matcher(text);
		StringBuffer sb = new StringBuffer();
		while (m.find())
		{
			Property property = Property.getById(Long.valueOf(m.group(1)));
			m.appendReplacement(sb, "<a tab=\"Property profile\" href=\"properties/edit.do?id="+property.id+"\">"+property.getName()+"</a>");
		}
		m.appendTail(sb);
		return sb.toString();
	}
	
	public static void main(String[] args) {
		System.out.println(transform(" ABC CND$1233 KJHB"));
	}
}
