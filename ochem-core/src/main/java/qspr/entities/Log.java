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

import java.text.Format;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlValue;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

@XmlRootElement(name = "log")
public class Log 
{
	@XmlTransient
	Format formatter = new SimpleDateFormat("HH:mm:ss");
	
	@XmlAttribute
	String name;
	
	@XmlAttribute
	Long id;
	
	@XmlAttribute
	String time;
	
	@XmlElement(name = "log-line")
	private List<LogLine> lines = new ArrayList<LogLine>();
	
	@XmlElement(name = "log")
	private List<Log> logs = new ArrayList<Log>();
	
	public Log addLog(String name)
	{		
		Log embeddedLog = new Log();
		embeddedLog.name = name;
		logs.add(embeddedLog);
		return embeddedLog;
	}
	
	
	public void addLine(String line, String type)
	{
		LogLine ll = new LogLine(line, formatter);
		ll.type = type;
		lines.add(ll);				
	}
	
	public void addLine(String line)
	{
		addLine(line, "normal");
	}
	
	public void addError(String line)
	{
		LogLine ll = new LogLine(line, formatter);
		ll.type = "error";
		lines.add(ll);
	}
	
	public void clear()
	{
		lines.clear();
	}
	
	public Log()
	{
		id = Math.round(Math.random()*100000);
		time = formatter.format(Calendar.getInstance().getTime());
	}
}

@XmlRootElement(name = "log-line")
class LogLine
{
	private static transient final Logger logger = LogManager.getLogger(LogLine.class);
	@XmlAttribute
	String time;
	
	@XmlAttribute
	String type = "normal";
	
	@XmlValue
	String value;
	
	public LogLine(String msg, Format formatter)
	{
		value = msg;
		time = formatter.format(Calendar.getInstance().getTime());
		logger.info(time+": "+msg);
	}
	
	public LogLine()
	{
	
	}
}
