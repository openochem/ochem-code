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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.poi.openxml4j.exceptions.InvalidFormatException;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;
import org.hibernate.Criteria;
import org.hibernate.Query;
import org.hibernate.criterion.Order;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Article;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Property;
import qspr.entities.Unit;
import qspr.export.ExportableMolecule;
import qspr.export.ExportableSetConfiguration;
import qspr.workflow.utils.QSPRConstants;
import qspr.util.AccessChecker;
import qspr.util.ExportThread;
import qspr.util.ValueDistribution;
import qspr.util.unitconversion.UnitConversion;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.MAILERConstants;
import com.eadmet.utils.OCHEMUtils;

@SuppressWarnings("unchecked")
public class BasketService 
{
	public List<Property> listUsedConditions(Basket basket)
	{
		Criteria conditionsCriteria = Globals.session().createCriteria(BasketEntry.class)
				.add(Restrictions.eq("basket", basket))
				.createAlias("ep", "expProperty")
				.createAlias("expProperty.conditions", "epConditions")
				.createCriteria("epConditions.values", "condVals")
				.setProjection(Projections.groupProperty("condVals.property"));
		conditionsCriteria.createAlias("condVals.property", "cond").add(Restrictions.ne("cond.type", 2));
		return conditionsCriteria.list();
	}

	public Basket clone(Basket b)
	{
		b.getPrivileges().requestView();
		Basket clone = new Basket();
		clone.session = Globals.userSession();
		clone.user = Globals.userSession().user;
		clone.name = OCHEMUtils.getFilteredBasketName("Copy of " + b.name);
		Globals.session().saveOrUpdate(clone);
		Globals.session().createSQLQuery("insert into BasketEntry(basket_id, exp_property_id, exclude) select "+clone.id+", exp_property_id, exclude from BasketEntry where basket_id=:basketId").setParameter("basketId", b.id).executeUpdate();
		return clone;
	}

	@SuppressWarnings("rawtypes")
	public Basket addRecordsFromAnotherBasket(Basket target, Basket source, String action)
	{
		String sql;
		if ("exclude".equals(action))
			sql = "update BasketEntry set exclude=1 where basket_id=:basketId and exp_property_id in (:list)";
		else if ("include".equals(action))
			sql = "update BasketEntry set exclude=0 where basket_id=:basketId and exp_property_id in (:list)";
		else if ("delete".equals(action))
			sql = "delete from BasketEntry where basket_id=:basketId and exp_property_id in (:list)";
		else
			sql = "insert ignore into BasketEntry(basket_id, exp_property_id) select "+target.id+", exp_property_id from BasketEntry source where source.basket_id=:sourceBasketId";

		Query query = Globals.session().createSQLQuery(sql);
		if (sql.contains(":sourceBasketId"))
			query.setParameter("sourceBasketId", source.id);
		if (sql.contains(":basketId"))
			query.setParameter("basketId", target.id);
		if (sql.contains(":list"))
		{
			List ids = Globals.session().createSQLQuery("select exp_property_id from BasketEntry where basket_id=" + source.id).list();
			query.setParameterList("list", ids);
		}

		query.executeUpdate();
		target.markModified();
		return target;
	}

	//Convert to universal parsers?
	public Basket addRecordsFromXls(Basket target, File f, String action) throws IOException, InvalidFormatException
	{
		String name = f.getName().toLowerCase();
		if (!name.contains(QSPRConstants.EXCEL) && !name.contains(QSPRConstants.EXCEL_NEW))
			throw new UserFriendlyException("Can not work with files which are not XLS");

		InputStream inp = new BufferedInputStream(new FileInputStream(f));
		Workbook wb = WorkbookFactory.create(inp);
		Sheet sheet = wb.getSheetAt(0);

		List<Long> ids = new ArrayList<Long>();

		Row header = sheet.getRow(0);
		for (int i=0; i<header.getLastCellNum(); i++)
		{
			String colName = header.getCell(i).getRichStringCellValue().toString();
			if(colName.equalsIgnoreCase("RECORDID"))
			{
				for (int r = 1; r <= sheet.getLastRowNum(); r++)
				{
					Row row = sheet.getRow(r);
					if (row == null)
						break;

					String id=row.getCell(i).getRichStringCellValue().toString();
					if(id.equalsIgnoreCase("RECORDID"))continue;

					if(id.equals(""))
						break;
					ids.add(Long.valueOf(id.replaceAll("R", "")));
				}
				break;
			}	    
		}

		inp.close();

		if (ids.size() > 0)
		{
			if ("exclude".equals(action))
				target.excludeOrIncludeEntries(ids, true);
			else if ("include".equals(action))
				target.excludeOrIncludeEntries(ids, false);
			else if ("delete".equals(action))
				target.removeEntries(ids);
			else
				target.addEntries(ids);
		}

		return target; 
	}

	public ExportThread getExportThread(final Long basketId, ExportableSetConfiguration conf, String format)
	{
		ExportThread eThread = new ExportThread(format, conf)
		{
			@Override
			public void generateData()
			{
				int firstEntry = 0;
				Basket basket;
				// Fetch entries in batches of size 1,000 and restart the transaction after each batch / Midnighter on Jun 6, 2011
				while (true)
				{
					basket = Basket.getBasket(Globals.userSession(), basketId);
					if(basket == null)throw new UserFriendlyException("This basket is not allowed to be accessed by OCHEM. If this is an error, contact us at "+MAILERConstants.EMAIL_OCHEM);
					basket.getPrivileges().requestView();
					List<BasketEntry> entries = Globals.session().createCriteria(BasketEntry.class).add(Restrictions.eq("basket", basket)).addOrder(Order.asc("id")).setMaxResults(1000).setFirstResult(firstEntry).list();
					firstEntry += entries.size();

					if (!entries.isEmpty())
					{
						for (BasketEntry entry : entries)
						{
							ExportableMolecule eMol = new ExportableMolecule();
							eData.addMolecule(eMol);
							eMol.setExperimentalProperty(entry.ep);
						}
						Globals.restartAllTransactions(true);
					}

					if (entries.size() < 1000)
						break;

					if (eData.exportableMolecules.size() % 100 == 0)
						setStatus("Preparing item " + eData.exportableMolecules.size());
				}
				setFileName(basket.name);
			}
		};
		return eThread;
	}

	public Basket getPrimaryBasket(Basket basket)
	{
		Basket primaryBasket = new Basket();
		primaryBasket.name = OCHEMUtils.getFilteredBasketName(basket.name + " (primary records)");
		primaryBasket.session = Globals.userSession();
		primaryBasket.user = Globals.userSession().user;
		Globals.session().saveOrUpdate(primaryBasket);
		Globals.session().flush();

		Set<Long> firstEntries = new HashSet<Long>();

		for(BasketEntry entry:basket.entries) 
		{
			if(entry.ep.firstEntry == null) {
				primaryBasket.entries.add(new BasketEntry(entry.ep)); // using same reference
				continue;
			}

			if(!firstEntries.contains(entry.ep.firstEntry)) {
				firstEntries.add(entry.ep.firstEntry);
				ExperimentalProperty ep = Repository.record.getRecord(entry.ep.firstEntry);
				if(AccessChecker.requestAddToBasket(ep))
					primaryBasket.entries.add(new BasketEntry(ep));
				else
					primaryBasket.entries.add(new BasketEntry(entry.ep));
			}
		}

		for(BasketEntry entry:primaryBasket.entries)
			entry.basket = primaryBasket;

		Globals.session().saveOrUpdate(primaryBasket);

		//Globals.session().createSQLQuery("insert ignore into BasketEntry(basket_id, exp_property_id) select " + primaryBasket.id + ", first_entry from BasketEntry natural left join ExperimentalProperty where basket_id=" + basket.id + " and first_entry is not null").executeUpdate();


		return primaryBasket;
	}

	public Basket fillDiscretizeMetadata(Basket basket)
	{
		String query = "select p.* from Basket b  left join BasketEntry be using (basket_id) left join ExperimentalProperty ep using (exp_property_id) left join Property p using (property_id) where b.basket_id = :id group by p.property_id";
		List<Property> properties = Globals.session().createSQLQuery(query).addEntity(Property.class).setLong("id", basket.id).list();

		basket.modelProperties = properties;
		for (Property property : properties) 
		{
			if (property.isNumeric())
			{
				List<Object[]> rows = Globals.session().createCriteria(BasketEntry.class)
						.createAlias("ep", "epa")
						.createAlias("epa.molecule", "mol")
						.add(Restrictions.eq("basket", basket))
						.add(Restrictions.eq("epa.property", property))
						.setProjection(Projections.projectionList().add(Projections.property("epa.value")).add(Projections.property("epa.unit")).add(Projections.property("mol.molWeight")))
						.list();

				List<Double> convertedValues = new ArrayList<Double>();
				for (Object[] row : rows)
					convertedValues.add(UnitConversion.convert((Double)row[0], (Unit)row[1], property.defaultUnit, (Double) row[2]));

				if (!convertedValues.isEmpty())
				{
					ValueDistribution dist = new ValueDistribution();
					dist.setValues(convertedValues);
					property.distribution = dist;
				}
			}
		}
		return basket;
	}

	public Basket fillEditMetadata(Basket basket)
	{
		basket.getPrivileges().requestView();
		basket.countRecordsByProperties(false); // count all the records
		basket.countRecordsByProperties(true); // count the excluded records
		basket.cachedCount = null;

		// Count the unique compounds
		Criteria uniqueCriteria = Globals.session().createCriteria(BasketEntry.class);
		uniqueCriteria.add(Restrictions.eq("basket", basket));
		uniqueCriteria.createAlias("ep", "e");
		uniqueCriteria.createAlias("e.molecule", "m");
		uniqueCriteria.createAlias("m.mapping1", "mp1");
		uniqueCriteria.createAlias("m.mapping2", "mp2");
		uniqueCriteria.setProjection(Projections.projectionList()
				.add(Projections.countDistinct("mp1.id"))
				.add(Projections.countDistinct("mp2.id")));

		List<Object[]> uniqueCompunds = uniqueCriteria.list();

		if(uniqueCompunds.size() > 0)
		{
			Object[] objects = uniqueCompunds.get(0);
			basket.totalUniqueCompounds = (Long) objects[0];
			basket.totalUniqueCompounds_SC = (Long)objects[1]; 
		}

		// count the records by articles
		Criteria articleCriteria = Globals.session().createCriteria(BasketEntry.class);
		articleCriteria.add(Restrictions.eq("basket", basket));
		articleCriteria.createAlias("ep", "e");
		articleCriteria.setProjection(Projections.projectionList().add(Projections.groupProperty("e.article")).add(Projections.count("id")));
		List<Object[]> artilceList = articleCriteria.list();
		if(artilceList.size() > 0)
		{
			basket.articleUsed = new ArrayList<Article>();
			for (Object[] objects : artilceList)
			{
				Article article = (Article) objects[0];
				article.count = (Long)objects[1];
				basket.articleUsed.add(article);
			}
		}

		// Count the records by tags
		Criteria tagCriteria = Globals.session().createCriteria(BasketEntry.class);
		tagCriteria.add(Restrictions.eq("basket", basket));
		tagCriteria.createAlias("ep", "e");
		tagCriteria.createAlias("e.molecule", "mol");
		tagCriteria.createAlias("mol.mapping1", "mp1");
		tagCriteria.createAlias("mp1.tags", "tag");
		tagCriteria.setProjection(Projections.projectionList().add(Projections.groupProperty("tag.id")).add(Projections.count("id")));
		/*
		List<Object[]> tagList = tagCriteria.list();
		if(tagList.size() > 0)
		{
			basket.tagUsed = new ArrayList<Tag>();
			for (Object[] objects : tagList)
			{
				Long id = (Long)objects[0];
				Tag tag = (Tag) Globals.session().get(Tag.class, id);
				try
				{
					AccessChecker.requestViewingPermission(tag);
					tag.count = (Long)objects[1]; 
					basket.tagUsed.add(tag);
				} catch (Exception e)
				{
					logger.info("Skipping non-public tag "+tag.name+" for molecules in basket "+basket.name);
				}
			}
		}
		 */
		return basket;
	}
	//private static Logger logger = Logger.getLogger(BasketService.class);
}
