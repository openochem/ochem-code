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

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlTransient;

public class XmlDate
{
	@XmlElement
	private int year;
	
	@XmlElement
	private int month;
	
	@XmlElement
	private int day;
	
	@XmlAttribute 
	public String DateType;
	
	@XmlTransient
	Date date;
	
	@XmlAttribute(name = "printed-name")
	public String getPrintedName()
	{
		if (date == null)
			return "";
		SimpleDateFormat df = new SimpleDateFormat("d MMM yyyy");
		return df.format(date);
	}
	
	public XmlDate() 
	{
	}
	
	public XmlDate(Date _date) 
	{
		date = _date;
		
		if (_date == null)
			return; 
		Calendar calendar = Calendar.getInstance();
		calendar.setTime(_date);
		
		year = calendar.get(Calendar.YEAR);
		month = calendar.get(Calendar.MONTH);
		day = calendar.get(Calendar.DAY_OF_MONTH);
	}

}
