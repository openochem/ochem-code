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

package qspr.metaserver.protocol;

import java.io.Serializable;
import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

/*
 * When server registers itself at Metaserver
 * It sends this class as a self-description
 * 
 * Midnighter
 */

public class ServerInfo implements Serializable
{
	private static final long serialVersionUID = 1L;
	public String platform;

	// Tasks supported by the server
	public List<String> supportedTaskTypes;

	// Tasks that failed to pass tests
	public Set<String> failures;

	public String name;
	public String owner;
	public String workingDirectory;
	public String version;
	public String status;
	public String ipAddress;

	/**
	 * Minimum task priority acceptable by this server
	 */
	public Integer minimumPriority;
	public long random;

	/**
	 * Available disc space in mbytes
	 */
	public long diskSpaceLeft;

	/**
	 * Currently used RAM in mbytes
	 */
	public long usedMemory;

	/**
	 * Available RAM in mbytes
	 */
	public long availableMemory;

	/**
	 * Approximate usage of CPU of the calculation server
	 */
	public double cpuUsage;

	/**
	 * Maxmimum RAM in mbytes
	 */
	public int peakMemory; 

	public String configurationXml;


	public void setSupportedTasks(List<String> tasks, String localConfig)
	{
		supportedTaskTypes = tasks;
		Collections.sort(supportedTaskTypes);
	}

	public void addFailure(String failure)
	{
		if (failures == null)
			failures = new TreeSet<String>();
		if (!failures.contains(failure))
			failures.add(failure);
	}

	public String toString()
	{
		return "" + failures + "\n" + status + "\n" + version;
	}
}
