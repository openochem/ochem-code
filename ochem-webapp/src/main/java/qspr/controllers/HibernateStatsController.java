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

package qspr.controllers;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.hibernate.stat.EntityStatistics;
import org.hibernate.stat.Statistics;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.entities.Log;
import qspr.frontend.WebModel;

@Controller
public class HibernateStatsController extends ControllerWrapper
{
	public ModelAndView list(HttpServletRequest request, HttpServletResponse response)
	{
		Log log = new Log();	
		Runtime.getRuntime().gc();
		Statistics stat = Globals.getHibernateStatistics();
	    log.addLine("Opened Sessions: "+stat.getSessionOpenCount());
	    log.addLine("Closed Sessions: "+stat.getSessionCloseCount());
	    log.addLine("Connections made :"+stat.getConnectCount());
	    log.addLine("Transactions :"+stat.getTransactionCount());
	    log.addLine("Succesfull transactions: "+stat.getSuccessfulTransactionCount());
	    log.addLine("Queries:"+stat.getQueryExecutionCount());
	    log.addLine("Max query time: "+stat.getQueryExecutionMaxTime());	
	    log.addLine("Entities fetched: "+stat.getEntityFetchCount());
	    log.addLine("Entities eventually fetched: "+stat.getEntityLoadCount());
	    log.addLine("Entities inserted: "+stat.getEntityInsertCount());
	    log.addLine("Entities updated: "+stat.getEntityUpdateCount());
	    log.addLine("Total runtime memory: "+Runtime.getRuntime().totalMemory());
	    log.addLine("Total free memory: "+Runtime.getRuntime().freeMemory());
	    log.addLine("Total max memory: "+Runtime.getRuntime().maxMemory());
	    String[] en = stat.getEntityNames();
	    StringBuffer sb = new StringBuffer();
	    for (int i=0; i<en.length; i++)
	    {
	    	EntityStatistics entityStats = stat.getEntityStatistics(en[i]);
	    	sb.append(en[i]);
	    	sb.append("(Fetched:");
	    	sb.append(entityStats.getFetchCount());
	    	sb.append(",Loaded:");
	    	sb.append(entityStats.getLoadCount());
    		sb.append("), ");
	    }
	    log.addLine("Entity names: "+sb.toString());
		return new WebModel(log).getModelAndView();
	}
}
