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

import java.sql.Timestamp;
import java.util.Calendar;
import java.util.GregorianCalendar;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import javax.xml.bind.annotation.XmlRootElement;
import javax.xml.namespace.QName;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.frontend.WebModel;

public class NewsFeedController
{
	private static transient final Logger logger = LogManager.getLogger(NewsFeedController.class);

	// A prototype for news feed
	@SuppressWarnings("unchecked")
	public ModelAndView show(HttpServletRequest request, HttpServletResponse response)
	{
		GregorianCalendar gCal = new GregorianCalendar();
		gCal.add(Calendar.MONTH, -1);
		Timestamp oneMonthAgo = new Timestamp(gCal.getTimeInMillis());

		List<Object[]> rows = Globals.session().createSQLQuery("select login, p.name, count(*) from Action natural left join Session s natural left join User inner join ExperimentalProperty ep on (ep.exp_property_id=primary_key) inner join Property p on (ep.property_id=p.property_id) where action_time>=:dateFrom and table_name='ExperimentalProperty' and comm='Item has been created' and ep.rights=2 group by login, p.property_id order by s.user_id")
				.setParameter("dateFrom", oneMonthAgo)
				.list();

		Object user = null;
		String msg = "";
		for (Object[] objects : rows)
		{
			if (user != null && !user.equals(objects[0]))
			{
				logger.info("User " + user + " has uploaded " + msg);
				msg = "";
			}
			msg += "" + objects[2] + " records for " + objects[1] + ", ";
			user = objects[0];
		}

		if (user != null)
			logger.info("User " + user + " has uploaded " + msg);

		return null;
	}

	public static void main(String[] args) throws JAXBException
	{
		Globals.jaxbContext = Globals.createJAXBContext();
		Marshaller marshaller = Globals.jaxbContext.createMarshaller();
		marshaller.marshal(new WebModel().setObject(new JAXBElement<Test>(new QName("array"), Test.class, new Test())), System.out);
		//		try
		//		{
		//			Globals.startAllTransactions();
		//			new NewsFeedController().show(null, null);
		//		}
		//		finally
		//		{
		//			Globals.commitAllTransactions();
		//		}

	}

	@XmlRootElement(name = "test")
	public static class Test
	{
		public int a = 1;
		public Long b = 2L;
	}
}


