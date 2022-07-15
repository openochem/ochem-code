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
import javax.persistence.Id;

@Entity
public class StatisticsLog 
{
	@Column(name = "log_id")
	@Id
	@GeneratedValue
	public Long id;
	
	@Column
	public Timestamp date;
	
	@Column(name = "connections_per_second")
	public Long connectionsPerSecond;
	
	@Column(name = "lost_connections_per_second")
	public Long lostConnectionsPerSecond;
	
	@Column(name = "assigned_tasks")
	public Long assignedTasks;
	
	@Column(name = "memory_used")
	public Long memoryUsed;
	
	@Column(name = "online_servers")
	public Long onlineServers;
	
	@Column(name = "new_tasks")
	public Integer newTasks;
	
	@Column(name = "completed_tasks")
	public Integer completedTasks;
	
	@Column
	public Integer errors;
	
	@Column(name = "tasks_in_queue")
	public Integer tasksInQueue;
}
