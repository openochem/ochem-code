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

import java.sql.Timestamp;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.Id;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.bind.annotation.XmlTransient;


@Entity
@XmlRootElement(name = "testresult")
public class PeriodicTestResult
{
	@Id
	@GeneratedValue
	@Column(name = "ptr_id")
	@XmlAttribute
	public Long id;
	
	@Column
	@XmlTransient
	public String name;
	
	@Column
	@XmlAttribute
	public boolean succeeded;
	
	@Column(name="detailed_status")
	@XmlElement
	public String detailedStatus;
	
	@Column(name="start_time")
	@XmlTransient
	public Timestamp startTime;

	@Column(name="finish_time")
	@XmlTransient
	public Timestamp finishTime;
	
	@Column(name="test_type")
	@XmlAttribute
	public String testType;
	
	@XmlAttribute
	public String getName()
	{
		return name.split("\\(")[0].replace("Test", ""); //.toUpperCase();
	}
	
	public void setName(String name)
	{
		this.name = name;
	}
	
	@XmlAttribute
	public String getTime()
	{
		return startTime.toString();
	}
	
	/**
	 * @return Duration of the test in seconds
	 */
	@XmlAttribute
	public long getDuration()
	{
		if (finishTime == null || startTime == null)
			return -1;
		// duration in secs
		return (finishTime.getTime() - startTime.getTime()) / 1000;
	}
	
	@Column(name="class_name")
	@XmlAttribute
	public String className;
	
	@Column(name="method_name")
	@XmlAttribute
	public String methodName;
}