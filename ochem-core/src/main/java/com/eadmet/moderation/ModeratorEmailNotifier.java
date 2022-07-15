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

package com.eadmet.moderation;

import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.entities.Property;
import qspr.entities.User;
import qspr.util.DynaWrap;
import qspr.util.WrapperThread;

import com.eadmet.datacube.DataTree;
import com.eadmet.datacube.Subkonto;
import com.eadmet.utils.MAILERConstants;
import com.eadmet.utils.mailer.Email;
import com.eadmet.utils.mailer.Mailer;

@SuppressWarnings("unchecked")
public class ModeratorEmailNotifier
{
	private static final Logger logger = LogManager.getLogger(ModeratorEmailNotifier.class);

	public static void notifyAllModerators() {
		List<User> moderators = getModerators();
		for (User moderator : moderators)
		{
			notifyByEmail(moderator);
			Globals.restartAllTransactions(true);
		}
	}

	public static List<User> getModerators() {
		return Globals.session().createCriteria(Property.class).add(Restrictions.isNotNull("moderator")).setProjection(Projections.groupProperty("moderator")).list();
	}

	public static void notifyByEmail(User moderator) {
		if (!(moderator.isExtended())) {
			logger.error("Current implementation of user does not support email addresses:", moderator.getClass().toString());
			// FIXME:  we should take alternative action here
			return;
		}
		
		ModeratorService service = new ModeratorService();
		DataReportRequest reportRequest = new DataReportRequest();
		reportRequest.groupingOrder = new ArrayList<Subkonto>();
		reportRequest.groupingOrder.add(new Subkonto(Property.class, "Property"));
		reportRequest.groupingOrder.add(new Subkonto(User.getCurrentClass(), "Introducer"));
		reportRequest.countPrivateRecords = false;
		reportRequest.countApprovedRecords = false;
		reportRequest.countRejectedRecords = false;
		reportRequest.countModels = false;
		reportRequest.moderator = moderator;
		
		DynaWrap wrappedModerator = moderator.dynaWrapped();

		AwaitingDataReport report = service.getAwaitingDataReport(reportRequest);

		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);

		pw.println("Dear " + moderator.getFullName() + ",\n");
		pw.println("Thank you for being a moderator at OCHEM.\nWe would like to kindly remind you that there is some data on OCHEM that is awaiting your approval.\n");
		pw.println("<b>The summary of the data awaiting your approval</b>:\n");

		pw.print("<table>");
		boolean anyData = printSummary(report.dataTree, pw, 0);
		pw.print("</table>");
		if (!anyData)
		{
			logger.info("No data to moderate for " + moderator);
			return;
		}

		pw.println("\nThe moderation corner is available at " +OCHEMConfiguration.rootHost+ "/editor/show.do?moderator=1.\nOnce you approve the data, you and the introducers of these data will get bonus points upon approval.\n");
		pw.println("Thank you for your collaboration and your contributions to the chemoinformatics community.\n\nBest regards,");
		pw.println("OCHEM Team");
		pw.flush();

		Mailer.postMailSafely(new Email(wrappedModerator.getString("email"), "Data awaiting your moderator's decision", sw.toString()).useHTML());

	}

	public static boolean printSummary(DataTree<DataMetrics> node, PrintWriter pw, int depth) {
		if (node.metrics.awaitingApproval > 0)
		{
			String title = depth == 0 ? "Total" : node.getTitle();
			pw.print("<tr>");
			pw.println("<td>" + prefix(depth, "&nbsp;") + title +  " </td><td align='right'>" + node.metrics.awaitingApproval + " records</td>");
			pw.print("</tr>");
			for (DataTree<DataMetrics> child : node.children)
				printSummary(child, pw, depth + 1);
			return true;
		}
		return false;
	}

	private static String prefix(int n, String charz) {
		if (n <= 0)
			return "";
		StringWriter sw = new StringWriter();
		for (int i = 0; i < n; i++)
			sw.append(charz);
		return sw.toString();
	}

	/*
	private static String addSpaces(String s, int n) {
		return s + prefix(n - s.length(), " ");
	}
	 */

	public static void main(String[] args)
	{
		new WrapperThread()
		{
			@Override
			public void wrapped() throws Exception
			{
				Mailer.enable = true;
				notifyByEmail(User.getByLogin(MAILERConstants.ADMIN));
			}
		}.run();
	}
}
