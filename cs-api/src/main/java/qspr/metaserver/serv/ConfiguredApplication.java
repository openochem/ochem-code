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

package qspr.metaserver.serv;

import java.util.Calendar;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlTransient;
import javax.xml.bind.annotation.XmlValue;

import qspr.metaserver.CalculationServer;

public class ConfiguredApplication 
{
	@XmlValue
	public String name;

	@XmlAttribute(name = "start")
	public Integer startTime;

	@XmlAttribute(name = "stop")
	public Integer stopTime;

	@XmlAttribute
	private Boolean disabled;

	@XmlTransient
	public CalculationServer server;

	public boolean isDisabled()
	{
		return disabled != null && disabled;
	}

	public void disable()
	{
		disabled = true;
	}

	public void enable()
	{
		disabled = null;
	}

	public boolean isAllowedToRun()
	{
		if (startTime == null || stopTime == null)
			return true;

		int currentHours = Calendar.getInstance().get(Calendar.HOUR_OF_DAY);
		if (startTime < stopTime)
			return currentHours >= startTime && currentHours <= stopTime;
			else
				return currentHours >= startTime || currentHours <= stopTime;
	}

	public ConfiguredApplication(String name)
	{
		this.name = name;
	}

	public ConfiguredApplication()
	{

	}
}
