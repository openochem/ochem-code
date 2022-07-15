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

package com.eadmet.utils;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.util.List;

import com.eadmet.exceptions.UserFriendlyException;

/**
 * Different utilities for gathering information about RAM usage
 * @author all
 */
public class MemoryUtils
{
	public final static float MEGABYTE = 1024*1024;

	/**
	 * Get the percentage of currently used memory (for the largest pool, which is usually OldGen) 
	 * @return current memory usage ratio
	 */
	public static double getCurrentMemoryUsedFraction()
	{
		MemoryPoolMXBean largestPool = getLargestPool();
		return Long.valueOf(largestPool.getUsage().getUsed()).doubleValue() / Long.valueOf(largestPool.getUsage().getMax()).doubleValue();
	}

	public static MemoryPoolMXBean getLargestPool()
	{
		MemoryPoolMXBean largestPool = null;
		List<MemoryPoolMXBean> pools = ManagementFactory.getMemoryPoolMXBeans();
		for(MemoryPoolMXBean pool : pools) 
			if (largestPool == null || pool.getUsage().getMax() > largestPool.getUsage().getMax())
				largestPool = pool;
		return largestPool;
	}

	public static void ensureSufficientMemory() {
		if (getCurrentMemoryUsedFraction() > 0.95)
			throw new UserFriendlyException("Unfortunately, your request cannot be processed because of insufficient memory on the OCHEM server. Please, try later or reduce the size of your dataset");
	}


	/**
	 * Get the current memory usage 
	 * @param oldGenerationOnly specifies whether we are interested only in the "Old generation" memory pool (usually, this is the case)
	 * @return current memory usage in MB
	 */
	public static int getCurrentMemoryUsage()
	{
		long result = 0;
		List<MemoryPoolMXBean> pools = ManagementFactory.getMemoryPoolMXBeans();
		for(MemoryPoolMXBean pool : pools) 
			result += pool.getUsage().getUsed();

		return Math.round(result / MEGABYTE);
	}


	/**
	 * Returns the textual summary of the memory pools usage 
	 * 
	 * @return
	 */
	public static String memorySummary()
	{
		StringBuffer sb = new StringBuffer();
		sb.append("[ ");
		List<MemoryPoolMXBean> pools = ManagementFactory.getMemoryPoolMXBeans();
		for(MemoryPoolMXBean pool : pools)
		{
			String[] pieces = pool.getName().split("\\s+");
			if (pieces.length > 1)
				sb.append(pieces[1].substring(0, 1));
			else
				sb.append("X");
		}
		sb.append("-UMA = ");
		for(MemoryPoolMXBean pool : pools) 
		{
			sb.append(NumericalValueStandardizer.getSignificantDigits(pool.getUsage().getUsed() / MEGABYTE));
			sb.append("mb. ");
			sb.append(NumericalValueStandardizer.getSignificantDigits(pool.getPeakUsage().getUsed() / MEGABYTE));
			sb.append("mb. ");
			sb.append(NumericalValueStandardizer.getSignificantDigits(pool.getUsage().getMax() / MEGABYTE));
			sb.append("mb + ");
		}
		sb.delete(sb.length()-3, sb.length());
		sb.append(" ]");
		return sb.toString();
	}

	/**
	 * Max heap VM can use e.g. Xmx setting
	 * @return
	 */

	public static int memoryVM(){
		Runtime runtime = Runtime.getRuntime();
		return Math.round(runtime.maxMemory() / MEGABYTE); // Max heap VM can use e.g. Xmx setting
	}

	/**
	 *  Available memory after deducing the one already occupied by Java
	 * @return
	 */

	public static int getCurrentMemoryFree()
	{
		int maxMemory = memoryVM(); // Max heap VM can use e.g. Xmx setting
		int usedMemory = getCurrentMemoryUsage(); // how much of the current heap the VM is using
		return maxMemory - usedMemory; // available memory i.e. Maximum heap size minus the current amount used
	}

	/**
	 * Get the peak memory usage 
	 * @param oldGenerationOnly specifies whether we are interested only in the "Old generation" memory pool (usually, this is the case)
	 * @return peak memory usage in MB
	 */
	public static int getPeakMemoryUsage()
	{
		long result = 0;
		List<MemoryPoolMXBean> pools = ManagementFactory.getMemoryPoolMXBeans();
		for(MemoryPoolMXBean pool : pools) 
			result += pool.getPeakUsage().getUsed();

		return Math.round(result / MEGABYTE); // how much of the current heap the VM is using
	}
	/*
	private static int getCurrentRuntimeMemoryUsage() {
		Runtime runtime = Runtime.getRuntime();
		long totalMemory = runtime.totalMemory(); // current heap allocated to the VM process
		long freeMemory = runtime.freeMemory(); // out of the current heap, how much is free
		return Math.round((totalMemory - freeMemory) / MEGABYTE); // how much of the current heap the VM is using
	}
	 */
}
