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

package qspr.metaserver;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/*
 * ServerPool
 * A singleton object containing available CalculationServers in given JVM
 * 
 *  Midnighter
 */

public class ServerPool 
{
	/**
	 * The only (singleton) instance of the Server Pool
	 */
	private static ServerPool instance = null;
	
	
	public String sid;
	public List<CalculationServer> servers = new ArrayList<CalculationServer>();
	
	public static ServerPool getInstance()
	{
		if (instance == null)
			instance = new ServerPool();
		return instance;
	}
	
	public CalculationServer getFreeServer(String taskType)
	{
		return getServer(taskType, true);
	}
	
	public CalculationServer getServer(String taskType, boolean freeServersOnly)
	{
		for (CalculationServer server : servers) {
			if (server.supportedTaskType.toLowerCase().equals(taskType.toLowerCase()) && (!server.busy || !freeServersOnly))
				return server;
		}
		
		return null;
	}
	
	public void disableTaskType(String taskType)
	{
		Iterator<CalculationServer> iterator = servers.iterator();
		while (iterator.hasNext())
			if (iterator.next().supportedTaskType.equals(taskType))
				iterator.remove();
	}

	
	/**
	 *  Check where calculations of this task can be performed locally
	 *  If yes, sets the task Id to negative value, which will be used to avoid (de)serilaisation
	 * @param task
	 * @return
	 */
	
	public static boolean canCalculateTaskLocally(String taskType)
	{
		return getInstance().getFreeServer(taskType) != null;
	}

}
