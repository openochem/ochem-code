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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.hibernate.Criteria;
import org.hibernate.FlushMode;
import org.hibernate.criterion.ProjectionList;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;
import org.springframework.stereotype.Service;

import qspr.Globals;
import qspr.entities.Basket;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.Property;
import qspr.entities.User;
import qspr.util.WrapperThread;

import com.eadmet.datacube.Subkonto;
import com.eadmet.datacube.SubkontoValue;

/**
 * A service for moderation-related activities
 * @author midnighter
 *
 */
@Service
@SuppressWarnings("unchecked")
public class ModeratorService 
{

	private List<Long> getRecordIds(Basket b, boolean approved)
	{
		Criteria c = Globals.session().createCriteria(ExperimentalProperty.class)
				.add(Restrictions.isNull("deleted"))
				.createAlias("basketEntries", "be")
				.add(Restrictions.eq("be.basket", b));

		if (!approved)
		{
			c.add(Restrictions.eq("approved", false));
			c.add(Restrictions.eq("rejected", false));
		} else
		{
			c.add(Restrictions.or(Restrictions.eq("approved", true), Restrictions.eq("rejected", true)));
		}

		c.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE))
		.setProjection(Projections.id());

		return c.list();
	}


	public void approveDataInBasket(Basket b)
	{
		List<Long> ids = getRecordIds(b, false);
		while (ids.size() > 0)
		{
			List<Long> subIds = ids.subList(0, Math.min(ids.size(), 100));
			Globals.session().createQuery("update ExperimentalProperty set approved = true where id in (:l)").setParameterList("l", subIds).executeUpdate();
			subIds.clear();
			Globals.restartAllTransactions(true);
		}

	}

	public void rejectDataInBasket(Basket b)
	{
		List<Long> ids = getRecordIds(b, false);
		while (ids.size() > 0)
		{
			List<Long> subIds = ids.subList(0, Math.min(ids.size(), 100));
			Globals.session().createQuery("update ExperimentalProperty set rejected = true, rights = 0 where id in (:l)").setParameterList("l", subIds).executeUpdate();
			subIds.clear();
			Globals.restartAllTransactions(true);
		}
	}

	public void unapproveDataInBasket(Basket b)
	{
		List<Long> ids = getRecordIds(b, true);
		while (ids.size() > 0)
		{
			List<Long> subIds = ids.subList(0, Math.min(ids.size(), 100));
			Globals.session().createQuery("update ExperimentalProperty set rejected = false, approved = false where id in (:l)").setParameterList("l", subIds).executeUpdate();
			subIds.clear();
			Globals.restartAllTransactions(true);
		}
	}

	/**
	 * Generate a report about the unapproved data
	 * 
	 * @param properties
	 * @param showApprovedRecords
	 * @return
	 */
	public AwaitingDataReport getAwaitingDataReport(DataReportRequest request)
	{
		AwaitingDataReport report = new AwaitingDataReport();
		report.cube.setupGroupingOrder(request.groupingOrder);
		List<Property> properties = null;

		if (request.moderator != null || !request.properties.isEmpty())
			properties = new ArrayList<Property>();

		if (request.moderator != null)
			properties.addAll(getModeratedProperties(request.moderator));
		if (!request.properties.isEmpty())
			properties.addAll(request.properties);

		Globals.session().setFlushMode(FlushMode.MANUAL);

		// Count unapproved records
		if (request.countUnapprovedRecords)
			addRecordData(report, properties, MetricsType.awaitingApproval, request);

		// Count approved records
		if (request.countApprovedRecords)
			addRecordData(report, properties, MetricsType.approvedRecords, request);

		if (request.countPrivateRecords)
		{
			addRecordData(report, properties, MetricsType.privateRecords, request);
			if (request.countModels)
				addModelsData(report, 0, request, properties);
		}

		if (request.countModels)
		{
			addModelsData(report, 1, request, properties);
			addModelsData(report, 2, request, properties);
		}

		if (request.moderator != null && request.countRejectedRecords)
			addRecordData(report, properties, MetricsType.rejected, request);


		Globals.session().setFlushMode(FlushMode.AUTO);



		report.dataTree = report.cube.getTree();

		return report;
	}

	private void addModelsData(AwaitingDataReport report, int type, DataReportRequest request, List<Property> properties) {
		// Awaiting models
		Criteria c = Globals.session().createCriteria(Model.class)
				.createAlias("modelMappings", "mm")
				.createAlias("session", "s")
				.setProjection(
						Projections.projectionList()
						.add(Projections.groupProperty("mm.property"))
						.add(Projections.groupProperty("s.user"))
						.add(Projections.countDistinct("id")));

		c.add(Restrictions.eq("published", type > 0));
		c.add(Restrictions.isNull("taskId"));
		if (type > 0)
			c.add(Restrictions.eq("approved", type == 2));	

		if (request.user != null)
			c.add(Restrictions.eq("s.user", request.user));

		if (properties != null)
			c.add(Restrictions.in("mm.property", properties));

		List<Object[]> rows = c.list();
		for (Object[] row : rows) 
		{
			DataMetrics metrics = new DataMetrics();
			int val = ((Long) row[2]).intValue();
			if (type == 0)
				metrics.privateModels = val;
			else if (type == 1)
				metrics.awaitingModels = val;
			else if (type == 2)
				metrics.approvedModels = val;
			report.cube.addValue(metrics, row[0], row[1]);
		}
	}

	private void addRecordData(AwaitingDataReport report, List<Property> properties, MetricsType metricsType, DataReportRequest request)
	{
		Criteria c = Globals.session().createCriteria(ExperimentalProperty.class).add(Restrictions.isNull("deleted"));

		if (request.qualitativeCondition != null)
		{
			c.createAlias("conditions", "cs", Criteria.LEFT_JOIN);
			c.createAlias("cs.values", "pv", Criteria.LEFT_JOIN, Restrictions.eq("pv.property", request.qualitativeCondition));
		}

		ProjectionList pList = Projections.projectionList();

		for (Subkonto subkonto : request.groupingOrder)
		{
			pList.add(Projections.groupProperty(subkonto.sqlGrouping));
		}
		pList.add(Projections.countDistinct("id"));

		c.setProjection(pList);


		if (request.user != null)
			c.add(Restrictions.eq("introducer", request.user));

		if (request.basketId != null)
		{
			Basket b = (Basket)Globals.session().get(Basket.class, request.basketId);
			c.createAlias("basketEntries", "be");
			c.add(Restrictions.eq("be.basket", b));
		}

		if (request.original != null)
		{
			if (request.original)
				c.add(Restrictions.eqProperty("id", "firstEntry"));
			else
				c.add(Restrictions.neProperty("id", "firstEntry"));
		}

		switch (metricsType) {
		case privateRecords:
			c.add(Restrictions.eq("rights", Globals.RIGHTS_NONE));
			break;
		case approvedRecords:
			c.add(Restrictions.eq("approved", true));
			c.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE));
			break;
		case awaitingApproval:
			c.add(Restrictions.eq("approved", false));
			c.add(Restrictions.eq("rejected", false));
			c.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE));
			break;
		case rejected:
			c.add(Restrictions.eq("approved", false));
			c.add(Restrictions.eq("rejected", true));
			c.add(Restrictions.eq("rights", Globals.RIGHTS_FREELY_AVAILABLE));
			break;
		case awaitingApprovalOverdue:
			break;
		case awaitingModels:
			break;
		default:
			break;
		}

		if (properties != null)
		{
			if (properties.isEmpty())
				return;
			c.add(Restrictions.in("property", properties));
		}


		// Awaiting records
		List<Object[]> rows = c.list();

		for (Object[] row : rows) {
			// No nulls supported yet
			//if (Arrays.asList(row).contains(null))
			//	continue;

			DataMetrics metrics = new DataMetrics();
			try
			{
				metrics.getClass().getField(metricsType.name()).set(metrics, ((Long) row[row.length - 1]).intValue());
			} catch (Exception e)
			{
				throw new RuntimeException(e);
			}

			Object[] subkontoValues = Arrays.copyOfRange(row, 0, row.length - 1);

			for (int i = 0; i < subkontoValues.length; i++)
				if (subkontoValues[i] == null)
					subkontoValues[i] = new SubkontoValue(request.groupingOrder.get(i), null);
			//for (int i = 0; i < subkontoValues.length)
			report.cube.addValue(metrics, subkontoValues);
		}
	}

	public List<Property> getModeratedProperties(User user)
	{
		return Globals.session().createCriteria(Property.class).add(Restrictions.eq("moderator", user)).list();
	}

	public static void main(String[] args)
	{
		new WrapperThread()
		{
			@Override
			public void wrapped() throws Exception
			{
				//AwaitingDataReport report = new ModeratorService().getAwaitingDataReport(null, true, new );
				//DataTree<DataMetrics> tree = report.cube.getTree();
				//Globals.createMarshaller(true).marshal(tree, System.out);
			}
		}.run();
	}
}

enum MetricsType {
	approvedRecords,
	privateRecords,
	awaitingApproval,
	awaitingApprovalOverdue,
	rejected,
	awaitingModels
}
