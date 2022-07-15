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

package qspr.modelling;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hibernate.Criteria;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.type.IntegerType;
import org.hibernate.type.Type;

import qspr.Globals;
import qspr.entities.PendingTask;
import qspr.entities.Session;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.protocol.Task.TaskPriority;

public class TaskPriorityManager 
{
	public static Map<Integer, Integer> getUserTotalQuota(Session s)
	{
		Map<Integer, Integer> map = new HashMap<Integer, Integer>();
		if (s == null || s.user == null) 
		{
			map.put(TaskPriority.EXTRA_HIGH, 0);
			map.put(TaskPriority.HIGH, 1);
			map.put(TaskPriority.NORMAL, 1);
			map.put(TaskPriority.LOW, 100);
		} else
			if (s.user.isOCHEMDeveloper())
			{
				map.put(TaskPriority.EXTRA_HIGH, 10);
				map.put(TaskPriority.HIGH, 100);
				map.put(TaskPriority.NORMAL, 1000);
				map.put(TaskPriority.LOW, 10000);
			} else
				if (s.user.isValidated())
				{
					map.put(TaskPriority.EXTRA_HIGH, 0);
					map.put(TaskPriority.HIGH, 10);
					map.put(TaskPriority.NORMAL, 100);
					map.put(TaskPriority.LOW, 1000);
				} else
				{
					map.put(TaskPriority.EXTRA_HIGH, 0);
					map.put(TaskPriority.HIGH, 1);
					map.put(TaskPriority.NORMAL, 10);
					map.put(TaskPriority.LOW, 100);
				}
		return map;
	}

	@SuppressWarnings("unchecked")
	private static Map<Integer, Integer> getUserFreeQuota(Session s, Long ignoredPendingTaskId)
	{
		// TODO
		// It was strange work-around since submission of tasks with conditions
		// required restarting transactions and resulted in exception here

		Map<Integer, Integer> map = getUserTotalQuota(s);
		try{
			Criteria c = Globals.session().createCriteria(PendingTask.class);
			if (s == null || s.user == null)
				c.add(Restrictions.eq("session", s));
			else
				c.createAlias("session", "s").add(Restrictions.eq("s.user", s.user));
			c.add(Restrictions.or(Restrictions.eq("status", Task.INIT), Restrictions.eq("status", Task.ASSIGNED), Restrictions.eq("status", Task.STOP)));
			c.add(Restrictions.isNotNull("priority"));
			if (ignoredPendingTaskId != null)
				c.add(Restrictions.ne("id", ignoredPendingTaskId));
			c.setProjection(Projections.sqlGroupProjection("priority p, count(*) c", "priority", new String[]{"p", "c"}, new Type[]{IntegerType.INSTANCE, IntegerType.INSTANCE}));
			List<Object[]> l = (List<Object[]>)c.list();
			for (Object[] res : l) 
			{
				Integer priority = (Integer)res[0];
				Integer count = (Integer)res[1];
				if (map.get(priority) != null)
					map.put(priority, map.get(priority) - count);
			}
		}catch(Exception ee){

		}
		return map;
	}

	public static Integer getNewTaskPriority(Integer requestedPriority, PendingTask p, Session s)
	{
		if(requestedPriority == TaskPriority.LOW)return TaskPriority.LOW;

		if (s != null && s.disableQuota)
			return requestedPriority;
		Map<Integer, Integer> map = getUserFreeQuota(s, p != null ? p.id : null);
		int currentPriority = requestedPriority;

		while (true)
		{
			if ((map.get(currentPriority) != null && map.get(currentPriority) > 0) || (currentPriority == TaskPriority.LOW))
				return currentPriority;
			currentPriority = TaskPriority.lowerPriority(currentPriority);
		}
	}

}
