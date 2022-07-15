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

import java.sql.Timestamp;

import javax.persistence.Column;
import javax.persistence.Entity;
import javax.persistence.GeneratedValue;
import javax.persistence.GenerationType;
import javax.persistence.Id;

import org.hibernate.annotations.GenericGenerator;


import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.TaskPriority;

//N.B.! Should be copy; cannot reuse the Archived due to Hibernate problems

@Entity
public class ArchivedTask {

	//N.B.! Should be copy; cannot reuse the Archived due to Hibernate problems


	@Column(name = "task_id", updatable = false, nullable = false)
	@GeneratedValue(strategy = GenerationType.AUTO, generator = "native")
	@GenericGenerator(name = "native", strategy = "native")
	@Id
	public Integer id;

	@Column(name = "reference_id")
	public String referenceId;

	@Column(name = "parent_task_id")
	public Integer parentTaskId;

	@Column
	public String status=null; // TODO Midnighter: put an enumerated type here

	@Column
	public String client;

	@Column
	public Integer datarows;

	@Column
	public Integer cols;

	@Column
	public Integer nonzero;

	@Column(name = "detailed_status")
	protected String detailedStatus;

	@Column(name="task_type")
	public String taskType;

	@Column(name = "calc_server_id")
	public String calcServerId;

	@Column
	public Timestamp time;

	@Column(name = "time_assigned")
	public Timestamp timeAssigned;

	@Column(name = "time_completed")
	public Timestamp timeCompleted;

	@Column
	public double priority = TaskPriority.NORMAL;

	@Column(name = "peak_memory_usage")
	public Integer peakMemoryUsage;

	@Column(name = "task_name")
	public String taskName;

	/**
	 * These columns contain data or references to data
	 */

	@Column
	protected byte[] configuration;

	@Column
	protected byte[] data;

	@Column
	protected byte[] result;


	@Column(name = "last_access")
	public Timestamp lastAccess;

	@Column(name = "scheduled_kill")
	public boolean scheduledKill;

	// In mbytes
	@Column(name = "min_required_memory")
	public Integer minRequiredMemory;

	// Number of times task was resubmitted
	@Column(name = "resubmitted")
	public Integer resubmitted;

	/**
	 * How many clients are interested in this task?
	 * 0 - the task in non-cachable
	 */
	@Column(name = "ref_count")
	public int referenceCount;

	/**
	 * MD5 key of the task input (data + configuration).
	 * Used for the task-level caching
	 */
	@Column(name = "task_md5")
	public String md5;

	@Column(name = "preferred_server")
	public String preferredServer;

	/**
	 * Indicates that this is a debug task and extra logging information should be returned by the calculation servers.
	 * 0 - no debug
	 * 1 - debug errors
	 * 2 - debug all tasks
	 */
	@Column
	public int debug = 0;

	/**
	 * The user who originated the calculation
	 */
	@Column
	protected String user;


	public ArchivedTask(Task task) {
		id = task.id;
		parentTaskId = task.parentTaskId;
		status = task.status;
		client = task.client;
		datarows = task.datarows;
		cols = task.cols;
		nonzero = task.nonzero;
		detailedStatus = task.getDetailedStatus();
		taskType = task.taskType;
		taskName = task.taskName;
		calcServerId = task.calcServerId;
		time = task.time;
		timeAssigned = task.timeAssigned;
		timeCompleted = task.timeCompleted;
		priority = task.priority;
		referenceId = task.referenceId;
		user = task.getUser();
		peakMemoryUsage = task.peakMemoryUsage;
		md5 = task.md5;
		minRequiredMemory = task.getMinRequiredMemory();
		resubmitted = task.resubmitted;
	}

	public ArchivedTask() {

	}

	public String getDetailedStatus() {
		return detailedStatus;
	}

	public String getUser() {
		return user;
	}
}
