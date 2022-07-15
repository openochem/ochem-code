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

package qspr.schedule;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.lang.management.OperatingSystemMXBean;
import java.util.ArrayList;
import java.util.List;

import com.eadmet.utils.MemoryUtils;

public class MonitoringTask extends OchemCronjobTask {

	/**
	 * The monitors contain values with 2 seconds intervals (configured in the Quartz timer dispatcher-servlet.xml)
	 */
	public static List<Long> memoryMonitor = new ArrayList<Long>();
	public static List<Double> loadMonitor = new ArrayList<Double>();
	
	public static OperatingSystemMXBean osBean = ManagementFactory.getOperatingSystemMXBean();
	public static MemoryMXBean memBean = ManagementFactory.getMemoryMXBean();
	
	@Override
	public void log(String st)
	{
		// Silent please
	}
	
 	public void executeTask() throws Exception
 	{
		memoryMonitor.add(MemoryUtils.getLargestPool().getUsage().getUsed() / (1024 * 1024));
		loadMonitor.add(1.0d * Math.round(osBean.getSystemLoadAverage() * 10) / 10);
		
		if (memoryMonitor.size() > 10000000)
			memoryMonitor.clear();
		
		if (loadMonitor.size() > 10000000)
			loadMonitor.clear();
	}
 	
	public static void main(String[] args)
	{
		new MonitoringTask().executeInternal();
	}
}
