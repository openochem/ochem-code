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

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/*  A convenience utility to write arrays in R files
 *  Midnighter 
 */

public class RWriter 
{
	public String variablePrefix = "";
	
	Map<String, StringWriter> map = new HashMap<String, StringWriter>();
	OutputStream out;
	List<String> lists = new ArrayList<String>();
	public long size = 0L;
	
	public RWriter(OutputStream outStream)
	{
		this.out = outStream;
	}
	
	public void doAppend(StringWriter writer, String piece)
	{
		size += piece.length();
		writer.append(piece);
	}
	
	public void doAppend(String varName, String piece)
	{
		doAppend(getBuffer(varName), piece);
	}
	
	public void addString(String varName, String value)
	{
		doAppend(getBuffer(varName), "\"" + value + "\",");
	}
	
	public void addValue(String varName, String value)
	{
		doAppend(getBuffer(varName), value + ",");
	}
	
	public void addValue(String varName, Double value)
	{
		if (value != null)
			doAppend(getBuffer(varName), "" + value + ",");
		else
			addNA(varName);
	}
	
	public void addValue(String varName, Integer value)
	{
		doAppend(getBuffer(varName), "" + value + ",");
	}
	
	public void addValue(String varName, Long value)
	{
		doAppend(getBuffer(varName), "" + value + ",");
	}
	
	public void addNA(String varName)
	{
		doAppend(getBuffer(varName), "NA,");
	}
	
	public void addValue(String varName, Object value)
	{
		if (value == null)
			addNA(varName);
		else if (value instanceof Long)
			addValue(varName, (Long) value);
		else if (value instanceof Integer)
			addValue(varName, (Integer) value);
		else if (value instanceof Double)
			addValue(varName, (Double) value);
		else
			addString(varName, value.toString());
	}
	
	public void writeAndClose() throws IOException
	{
		write();
		close();
	}
	
	public void write()
	{
		PrintWriter pw = new PrintWriter(out);
		createLists(pw);
		for (String var : map.keySet())
		{
			pw.write(var + " = c(");
			pw.write(map.get(var).toString());
			pw.write("NULL)\n");
		}
		pw.flush();
		map.clear();
	}
	
	public void writeDirect(String line) throws IOException
	{
		out.write(("\n" + line + "\n").getBytes());
	}
	
	public void close() throws IOException
	{
		out.flush();
		out.close();
	}
	
	private StringWriter getBuffer(String name)
	{
		name = variablePrefix + name;
		name = name.replaceAll("\\s", ".");
		StringWriter writer = map.get(name);
		if (writer == null)
			map.put(name, writer = new StringWriter());
		return writer;
	}
	
	private void createLists(PrintWriter pw)
	{
		// Create all necessary lists
		for (String var : map.keySet())
		{
			if (var.contains("$"))
			{
				String[] parts = var.split("\\$");
				String list = "";
				for (int i = 0; i < parts.length - 1; i++)
				{
					list += parts[i];
					if (!lists.contains(list))
					{
						pw.write(list + " = list()\n");
						lists.add(list);
					}
					list += "$";
				}
			}
		}
	}
}
