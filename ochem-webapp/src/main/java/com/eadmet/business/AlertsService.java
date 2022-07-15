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

import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.sql.Timestamp;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.hibernate.HibernateException;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.business.toxalert.AlertsFilter;
import qspr.business.toxalert.ScreeningProcessor;
import qspr.dao.Repository;
import qspr.entities.Article;
import qspr.entities.ConditionSet;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.PropertyValue;
import qspr.entities.SubstructureAlert;
import qspr.export.CSVExportWriter;
import qspr.export.ExportableMolecule;
import qspr.export.ExportableSetConfiguration;
import qspr.metaserver.util.ExtendedSMART;
import qspr.util.ExportThread;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.parsers.SimpleParser;

public class AlertsService
{

	public static ExportThread getResultExporter(String format, ExportableSetConfiguration config, final ScreeningProcessor processor)
	{
		ExportThread eThread = new ExportThread(format, config)
		{
			@Override
			@Deprecated
			public void generateData() throws Exception
			{
				DataTable dtAlerts = new DataTable(true);

				//TODO delete old code after sometimes
				if(processor.orderedCompounds == null) // old code, to be deleted after six months, 17 January 2020
					for (int i = 0; i < processor.alertsByCompounds.index.size(); i++)
					{
						ExportableMolecule eMol = new ExportableMolecule();
						eData.addMolecule(eMol);
						eMol.setMolecule(Repository.molecule.getMapping2(processor.alertsByCompounds.index.get(i)).getMolecule());

						dtAlerts.addRow();
						eMol.descriptors = dtAlerts.getCurrentRow();
						for (int k = 0; k < processor.alertsByCompounds.values.get(i).size(); k++)
							dtAlerts.setValue("TA" + processor.alertsByCompounds.values.get(i).get(k), "1");

						if (i % 100  == 0)
							Globals.restartAllTransactions(true);

						setStatus("Prepared " + i + " records out of " + processor.alertsByCompounds.index.size());
					}
				else
				{ // new code
					int all = 0;

					for(Integer mol:processor.orderedCompounds) {

						ExportableMolecule eMol = new ExportableMolecule();
						eData.addMolecule(eMol);
						dtAlerts.addRow();
						if(mol == null) {
							eMol.error = "calculation failed";
							dtAlerts.getCurrentRow().setError(eMol.error);
							continue; // molecule failed
						}

						eMol.setMolecule(Repository.molecule.getMapping2(mol).getMolecule());

						Integer index = processor.alertsByCompounds.indexPositionHash.get(mol); 
						if(index == null) 
							continue; // no results were found, thus all values are zero

						eMol.descriptors = dtAlerts.getCurrentRow();
						for (int k = 0; k < processor.alertsByCompounds.values.get(index).size(); k++)
							dtAlerts.setValue("TA" + processor.alertsByCompounds.values.get(index).get(k), "1");

						if (all % 100  == 0)
							Globals.restartAllTransactions(true);

						setStatus("Prepared " + all + " records out of " + processor.alertsByCompounds.index.size());
					}
				}

				eData.setDescriptors(dtAlerts);
				setFileName("screening-results");
			}
		};
		return eThread;
	}

	public static AlertsUploadReport uploadAlerts(InputStream is, String fileName, boolean privateUpload, Property defaultEndpoint, Article defaultArticle) throws NumberFormatException, Exception {

		AlertsUploadReport report = new AlertsUploadReport();
		SimpleParser parser = SimpleParser.getParser(fileName).setSource(is);
		parser.reset();

		Map<SASheetColumn, Integer> colPosition = new HashMap<SASheetColumn, Integer>();

		Map<Property, Integer> conditionColums = new HashMap<Property, Integer>();

		List<String> headerCols = parser.sheetColumns.get(0);
		for (int i = 0; i < headerCols.size(); i++) {
			if (headerCols.get(i) == null)
				continue;
			String col = headerCols.get(i);
			try {
				colPosition.put(SASheetColumn.valueOf(col.toUpperCase()), i);
			} catch (Exception e) {
				Property prop = Repository.property.getProperty(col, false);
				if (prop == null)
					report.warnings.add("Unknown column " + col);
				else
					conditionColums.put(prop, i);
			}
		}

		while (parser.hasNext())
		{
			List<String> values = parser.next();
			try {
				SubstructureAlert alert = new SubstructureAlert();
				alert.rights = privateUpload ?  Globals.RIGHTS_NONE : Globals.RIGHTS_FREELY_AVAILABLE;
				alert.smart = values.get(colPosition
						.get(SASheetColumn.SMARTS));

				String property = values.get(colPosition.get(SASheetColumn.PROPERTY));
				if (property != null)
				{
					alert.property = Repository.property.getProperty(property, false);
					if (alert.property == null)
						throw new UserFriendlyException("Unknown property "
								+ values.get(colPosition
										.get(SASheetColumn.PROPERTY)) + " check property name by exact match to existing ones or create a new one before uploading the alerts");
				}
				else
					alert.property = defaultEndpoint;

				if ("".equals(alert.smart))
					continue;

				String articleStr = null;
				if (colPosition.get(SASheetColumn.PUBMEDID) != null)
				{
					articleStr = values.get(colPosition.get(SASheetColumn.PUBMEDID));
					if (articleStr != null)
					{
						String pmId = values.get(colPosition.get(SASheetColumn.PUBMEDID));
						if (!"".equals(pmId))
							alert.article = Article.getByPmId(new Long(pmId));
					}
				}

				if (colPosition.get(SASheetColumn.ARTICLEID) != null)
				{
					articleStr = values.get(colPosition.get(SASheetColumn.ARTICLEID));
					if (alert.article == null && articleStr != null)
						alert.article = Article.getArticle(articleStr);
				}

				if (alert.article == null)
					alert.article = defaultArticle;

				if (alert.article == null)
					throw new UserFriendlyException("Unknown article "
							+ values.get(colPosition
									.get(SASheetColumn.ARTICLEID)));
				else
					Globals.session().saveOrUpdate(alert.article);

				alert.introducer = alert.owner = Globals.userSession().user;
				alert.timeCreated = alert.time = new Timestamp(Calendar
						.getInstance().getTimeInMillis());

				ExtendedSMART exSmart = ExtendedSMART.create(alert.getFullSMARTS(), OCHEMConfiguration.getCheminfEngine());
				if (exSmart.invalid)
					throw new UserFriendlyException("Invalid SMART: "
							+ alert.smart);

				if (colPosition.containsKey(SASheetColumn.COMMENT))
					alert.comment = values.get(colPosition
							.get(SASheetColumn.COMMENT));
				if (colPosition.containsKey(SASheetColumn.NAME))
					alert.name = values.get(colPosition
							.get(SASheetColumn.NAME));
				if (colPosition.containsKey(SASheetColumn.SMARTS_DESCRIPTION))
					alert.smartsDescription = values.get(colPosition
							.get(SASheetColumn.SMARTS_DESCRIPTION));
				if (colPosition.containsKey(SASheetColumn.DESCRIPTION))
					alert.description = values.get(colPosition
							.get(SASheetColumn.DESCRIPTION));

				ConditionSet cSet = new ConditionSet();
				for (Entry<Property, Integer> condColumn : conditionColums
						.entrySet()) {
					String val = values.get(condColumn
							.getValue());
					if (val != null && !"".equals(val)) {
						PropertyOption option = condColumn.getKey()
								.getOptionByName(val);
						if (option == null)
							throw new UserFriendlyException("Unknown option "
									+ val + " for condition "
									+ condColumn.getKey().getName());
						PropertyValue pVal = new PropertyValue(condColumn
								.getKey().getOptionByName(val));
						cSet.values.add(pVal);
					}
				}
				alert.approved = false;

				if (!cSet.values.isEmpty())
					alert.conditions = cSet.get();

				alert.updateHash();

				if (!alert.hasConflicts()) {
					Globals.session().save(alert);
					report.successes++;
				} else
					report.dublicates++;
			} catch (Exception e) {
				e.printStackTrace();
				report.warnings.add(e.getMessage());
				report.errors++;
			}
		}

		return report;
	}

	@SuppressWarnings("unchecked")
	public static void exportAlerts(OutputStream os) throws HibernateException, Exception
	{
		SimpleDateFormat year = new SimpleDateFormat("yyyy");
		CSVExportWriter eWriter = new CSVExportWriter();
		eWriter.os = os;
		eWriter.pw = new PrintWriter(os);
		eWriter.initialize();
		eWriter.writeRow(new String[]{
				"Alert ID",
				SASheetColumn.SMARTS.toString(),
				SASheetColumn.ARTICLEID.toString(),  
				SASheetColumn.PROPERTY.toString(),
				SASheetColumn.NAME.toString(),
				SASheetColumn.DESCRIPTION.toString(),
				SASheetColumn.PAGE.toString(), 
				SASheetColumn.TABLE.toString(),
				SASheetColumn.LINE.toString(),
				SASheetColumn.COMMENT.toString(),
				SASheetColumn.SMARTS_DESCRIPTION.toString(),
				SASheetColumn.PUBMEDID.toString(),
		"Publication"});


		List<SubstructureAlert> alerts = new AlertsFilter().filterCriteria().list();
		for (SubstructureAlert alert : alerts)
			eWriter.writeRow(new Object[]{"TA" + alert.id, alert.smart, "A" + alert.article.id, alert.property.getName(), alert.name, alert.description, alert.artPageNum, alert.artTableNum, alert.artLineNum, alert.comment, alert.smartsDescription,
					alert.article.pmid, year.format(alert.article.publicationDate) + " " + alert.article.authors.get(0).lastName});

		eWriter.flush();
	}

	enum SASheetColumn {
		SMARTS, DESCRIPTION, SMARTS_DESCRIPTION, ARTICLEID, PROPERTY, NAME, COMMENT, PAGE, TABLE, LINE, PUBMEDID
	}
}
