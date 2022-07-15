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

package qspr.workflow.datatypes;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.protocol.Task;
import qspr.workflow.utils.QSPRConstants;

public class MemoryEstimator
{
	private static transient final Logger logger = LogManager.getLogger(MemoryEstimator.class);

	/**
	 * An estimate for the RAM required to complete this task.
	 * In Megabytes.
	 */
	public static int getMinRequiredMemory(Task t, WorkflowNodeData wnd)
	{
		if(t.isLocalTask()) return 0; // all local task are assumed to have sufficient memory by design

		if(wnd.getRows() < 50) return  768;  // few molecule only -- should be processed as fast as possible and should fit really small memory requirements
		if(wnd.getRows() <= 1077) return  1024;  // ALL our test tasks should fit to 1 GB ...

		// allocation based on the size of data
		int memoryDataSize = (int) Math.ceil(wnd.getRows()/30);  // for each 30 records we allocate 1MB of memory, e.g. for 30k records we allocate 1GB before the multiplication

		// allocation based on the size of model
		int memoryCfgSize = (int)Math.ceil(t.getConfigurationSize() / (128 * 1024)); // for each 128 KB we allocate 1MB of memory

		int memory = Math.max(memoryDataSize, memoryCfgSize) * memoryMultiplier(t.taskType, t.taskName); // allocates more memory for some specific tasks

		memory = memory < 512 ? 1024 : 1024 * (int)Math.ceil(0.5 + memory/1024.); // rounding to full MB

		switch(t.taskType) {
		case QSPRConstants.DLCA: memory = memory < 3072? 3072: memory;
		case QSPRConstants.CORINA:
		case QSPRConstants.DNN: 
		case QSPRConstants.PYTORCH: 
		case QSPRConstants.DEEPCHEM: 
		case QSPRConstants.LSSVMG: 
		case QSPRConstants.TRANSNN: 
		case QSPRConstants.CHEMPROP: 
		case QSPRConstants.CNF: memory = memory > 8192 ? 8192: memory; break;
		default: memory = memory > 65536 ? 65536: memory; break; // no more than 64GB for the moment
		}

		if(t.taskName != null && t.taskName.toUpperCase().contains(QSPRConstants.DLCA))
			memory = memory < 3072? 3072: memory;

		logger.info("Memory has been determined as "+memory +"MB (task type "+t.taskType+";" + t.taskName + ")");

		return memory;
	}

	/**
	 *  Allocates even more memory for very large tasks
	 */

	private static int memoryMultiplier(String taskType, String taskName){

		taskType = taskType == null ?"" : taskType.toLowerCase();
		taskName = taskName == null ?"" : taskName.toLowerCase();

		String workflowlarge[]= {QSPRConstants.Workflow,QSPRConstants.CONSENSUS};
		String workflowmedium[]= {QSPRConstants.DESCRIPTORS,QSPRConstants.BAGGING,"chemaxondescriptors",
				QSPRConstants.SELECTION,QSPRConstants.CROSSVALIDATION,QSPRConstants.PLS,QSPRConstants.J48};

		String largeMemories[] = {"pydescriptor", "dragon", "deepchem", "rdkit","sirms","macau","mmpfrag","chemaxonscaffold", "transnn", "epa"};

		for(String task:largeMemories) // this task itself
			if(task.equalsIgnoreCase(taskType)) return 5;

		for(String task:workflowlarge) // is part of workflow
			if(task.equalsIgnoreCase(taskType)) {
				for(String name:largeMemories)
					if(taskName.contains(name)) return 9; // checking whether we have very large memory requirements
				return 3;
			}

		for(String task:workflowmedium)
			if(task.equalsIgnoreCase(taskType)){return 2;}

		return 1;
	}


}
