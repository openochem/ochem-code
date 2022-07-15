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

package qspr.metaserver.frontend;


import java.util.List;

import org.hibernate.Transaction;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Restrictions;

import qspr.metaserver.ArchivedTask;
import qspr.metaserver.MetaServer;
import qspr.metaserver.protocol.Task;

public class QuickReport 
{
	public void printReport()
	{
		MetaServer ms = new MetaServer(MetaServer.METASERVER);
		System.out.println("Starting query");
		Transaction tx = ms.session().beginTransaction();
		@SuppressWarnings("unchecked")
		List<ArchivedTask> atl = 
		ms.session().createCriteria(ArchivedTask.class)
		.add(Restrictions.eq("taskType",Task.Workflow))
		.add(Restrictions.ge("rows", 10000))
		.add(Restrictions.eq("status", "ready"))
		.add(Restrictions.not(Restrictions.like("taskName", "Applying%")))
		.add(Restrictions.not(Restrictions.like("taskName", "%KNN%")))
		.add(Restrictions.like("taskName", "%Bag%"))
		.add(Restrictions.eq("user", "enamine"))
		//				.add(Restrictions.eq("id", 29287986))
		.addOrder(Order.desc("id"))
		.setMaxResults(1000)
		.list();
		for (ArchivedTask archivedTask : atl) 
			System.out.println(archivedTask.taskName);
		System.out.println("Done query, starting aggregation");
		for (ArchivedTask archivedTask : atl) 
		{
			TaskTreeNode node = new TaskTreeNode(archivedTask);
			String taskName = archivedTask.taskName;
			String validation = "CV";
			if (taskName.contains("Bag"))
				if (taskName.contains("SBag"))
					validation = "SBAG";
				else
					validation = "BAG";

			String[] pieces = taskName.split("_");
			if (pieces.length < 3)
				continue;
			System.out.println(String.format("%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s", archivedTask.id, archivedTask.datarows, node.maxColumns, MetaServer.niceTime(node.task.getCalculationTime()), MetaServer.niceTime(node.totalTime), pieces[0], validation, pieces[1], pieces[2]));
		}
		tx.rollback();
	}

	public static void main(String[] args)
	{
		new QuickReport().printReport();
	}
}
