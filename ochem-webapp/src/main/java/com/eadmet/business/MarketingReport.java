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

package com.eadmet.business;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.List;

import javax.xml.bind.annotation.XmlRootElement;

import org.hibernate.Criteria;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.hibernate.type.LongType;

import qspr.Globals;
import qspr.entities.User;
import qspr.util.WrapperThread;


@XmlRootElement
public class MarketingReport {
	public List<WeeklySummary> weeklySummary = new ArrayList<WeeklySummary>();

	public static MarketingReport create(int minActions) {
		MarketingReport report = new MarketingReport();

		for (int i = 0; i < 20; i++)
		{
			WeeklySummary summary = new WeeklySummary();

			Calendar cal = Calendar.getInstance();
			cal.setFirstDayOfWeek(Calendar.MONDAY);
			cal.add(Calendar.WEEK_OF_YEAR, i-20);

			summary.weekNum = cal.get(Calendar.WEEK_OF_YEAR);
			summary.yearNum = cal.get(Calendar.YEAR);
			List<User> newUsers = getRegisteredUsersByWeek(summary.weekNum, summary.yearNum, minActions);
			summary.newRegistrations = newUsers.size();
			for (User user : newUsers)
				summary.newRegistrationsLogins.add(user.login);

			report.weeklySummary.add(summary);
		}

		return report;
	}

	public static class WeeklySummary {
		public int weekNum;
		public int yearNum;
		public int newRegistrations;
		public List<String> newRegistrationsLogins = new ArrayList<String>();
	}

	@SuppressWarnings("unchecked")
	public static List<User> getRegisteredUsersByWeek(int week, int year, int minActions) {
		Calendar cal = Calendar.getInstance();
		cal.setFirstDayOfWeek(Calendar.MONDAY);
		cal.set(Calendar.YEAR, year);
		cal.set(Calendar.WEEK_OF_YEAR, week);
		cal.set(Calendar.DAY_OF_WEEK, Calendar.MONDAY);

		Date from = cal.getTime();
		cal.add(Calendar.WEEK_OF_YEAR, 1);
		Date to = cal.getTime();

		System.out.println(from);
		System.out.println(to);

		Criteria c = Globals.session().createCriteria(User.getCurrentClass(), "u").add(Restrictions.between("registrationTime", from, to));

		if (minActions > 0)
		{
			List<Long> userIDs = Globals.session().createSQLQuery("select User.user_id from User natural join Session natural join UserEvent group by user_id having count(*) > :minActions")
					.addScalar("user_id", LongType.INSTANCE)
					.setInteger("minActions", minActions)
					.list();
			c.add(Restrictions.in("id", userIDs));
		}

		c.setProjection(Projections.distinct(Projections.id()));

		List<Long> ids = c.list();

		if (ids != null && !ids.isEmpty())
			return Globals.session().createCriteria(User.getCurrentClass()).add(Restrictions.in("id", ids)).list();
		else
			return new ArrayList<User>();
	}

	public static void main(String[] args) {
		new WrapperThread() {

			@Override
			public void wrapped() throws Exception {
				System.out.println(MarketingReport.getRegisteredUsersByWeek(35, 2013, 5));
			}
		}.run();;
	}
}
