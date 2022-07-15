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

package qspr.metaserver.cs.util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.jezhumble.javasysmon.JavaSysMon;
import com.jezhumble.javasysmon.OsProcess;

public class CpuMeasurementTool 
{
	private static transient final Logger logger = LogManager.getLogger(CpuMeasurementTool.class);

	private static long MEASUREMENT_BUFFER = 10 * 60 * 1000;
	private static long MEASUREMENT_PRECISION = 1000;
	private static int[] MEASUREMENT_INTERVALS = {5 * 60 * 1000, 1 * 60 * 1000, 5 * 1000};
	class Measure
	{
		long processMillis = 0;
		long cpuMillis = 0;
		public Measure()
		{
			cpuMillis = Calendar.getInstance().getTimeInMillis();
			processMillis = 0;
			for (OsProcess proc : getProcessTreeList(sysmon.processTree().find(sysmon.currentPid()))) 
				processMillis += proc.processInfo().getSystemMillis()+proc.processInfo().getUserMillis();
		}

		public double getLoad(Measure m)
		{
			if (m.cpuMillis == cpuMillis)
				return 0.0;

			if (m.processMillis < processMillis)
				return 0.0;

			return (m.processMillis - processMillis) * 1.0 / (m.cpuMillis - cpuMillis);
		}
	}

	JavaSysMon sysmon = new JavaSysMon();
	List<Measure> measures = new ArrayList<Measure>();


	@SuppressWarnings("unchecked")
	private List<OsProcess> getProcessTreeList(OsProcess proc)
	{
		List<OsProcess> treeList = new ArrayList<OsProcess>();
		treeList.add(proc);

		int index = 0;
		while (index < treeList.size())
		{
			treeList.addAll(treeList.get(index).children());
			index++;
		}

		return treeList;
	}


	public Measure measure()
	{
		Measure m = new Measure();

		if (measures.size() == 0 || (m.cpuMillis - measures.get(measures.size()-1).cpuMillis > MEASUREMENT_PRECISION))
			measures.add(m);

		while (m.cpuMillis - measures.get(0).cpuMillis > MEASUREMENT_BUFFER)
			measures.remove(0);

		return measures.get(measures.size()-1);
	}


	public double[] getLoads()
	{
		double[] result = new double[MEASUREMENT_INTERVALS.length];
		Arrays.fill(result, 0);

		try {

			Measure m = measure();

			for (int index=0; index<MEASUREMENT_INTERVALS.length; index++)
			{
				Measure lastMeasure = null;
				for (Measure cm : measures) 
				{
					if (m.cpuMillis - cm.cpuMillis < MEASUREMENT_INTERVALS[index])
					{
						if (!m.equals(cm) || lastMeasure == null)
							result[index] = cm.getLoad(m);
						else
							result[index] = lastMeasure.getLoad(m);
						break;
					}
					lastMeasure = cm;
				}
			}
		} catch (Exception e)
		{
			logger.info("Exception in getLoads, possibly an attempted call on a killed process");
		}

		return result;
	}

	public static void main(String[] args)
	{
		double[] loads = new CpuMeasurementTool().getLoads();
		System.out.println("" + loads[0] + " " + loads[1]);
	}
}
