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

package qspr.util;

import java.sql.Timestamp;
import java.util.Calendar;

import qspr.Globals;
import qspr.entities.ExportAction;
import qspr.export.ExportWriter;
import qspr.export.ExportableSet;
import qspr.export.ExportableSetConfiguration;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.FileUtils;

/**
 * An abstract thread that prepares the exportable data
 * Based on the Long Operation API
 * 
 * @author midnighter/novserj
 *
 */
public abstract class ExportThread extends WrapperThread
{
	private String fileName;
	public long operationID;

	/**
	 * The directory for the generated file
	 */
	public ExportAction eAction;

	/**
	 * The dataset
	 */
	protected ExportableSet eData;

	String format;

	public ExportThread(String format, ExportableSetConfiguration conf)
	{
		eData = new ExportableSet();
		operationID = Operation.generateID();
		this.format = format;
		eData.configure(conf);
	}

	/**
	 * The method that actually prepares the dataset (problem-specific)
	 * @throws Exception
	 */
	public abstract void generateData() throws Exception;


	public void setFileName(String name)
	{
		String replacePattern = "[\\s;\\:\\-()=/<>.,\\\\+_]+";
		fileName = name.replaceAll(replacePattern, "_");
	}

	public void wrapped() throws Exception
	{
		registerExportEvent();
		Globals.restartAllTransactions(true);
		setOperationID(operationID);


		try
		{
			generateData();

			// Automatically replace Excel with CSV if there are more rows than Excel allows
			if (QSPRConstants.EXCEL.equals(format) && eData.exportableMolecules.size() >= 65534)
				eData.exportableMolecules = eData.exportableMolecules.subList(0, 65534);
			//format = "csv";

			eAction.setFileName(fileName);
			fileName = eAction.getFileName(); //The eAction cuts the filename length to 98 characters... we need to take it into account
			fileName = FileUtils.convertToASCII(fileName);
			eAction.setFileName(fileName);
			Globals.session().saveOrUpdate(eAction);
			//Update the overrides that could have happened in the meantime
			Globals.restartAllTransactions(true);

			setStatus("Writing items...");
			ExportWriter eWriter = ExportWriter.createWriter(format, eData, eAction.getExportFolder(), fileName);
			eAction.format = eWriter.getFileExtension();
			Globals.session().saveOrUpdate(eAction);

			operation.successURL = "/export/download.do?id=" + eAction.id;
			eWriter.write();
			setStatus("Finished");
			finalizeExportEvent(eWriter);
		} catch (Exception e)
		{
			eAction.failed = true;
			eAction.status = e.getMessage();
			Globals.session().saveOrUpdate(eAction);
			Globals.restartAllTransactions(true);
			throw e;
		}
	}

	private void registerExportEvent()
	{
		eAction = new ExportAction();
		eAction.session = Globals.userSession();
		eAction.exportedColumns = eData.selectedColumns.toString();
		eAction.format = format;
		eAction.status = "Initializing";
		eAction.time = new Timestamp(Calendar.getInstance().getTimeInMillis());
		Globals.session().saveOrUpdate(eAction);
	}

	private void finalizeExportEvent(ExportWriter eWriter)
	{
		eAction.count = eData.exportableMolecules.size();
		eAction.totalBonusPoints = eWriter.exportCost;
		if (eAction.totalBonusPoints == 0) //Free downloads are automatically accepted
			eAction.fileDownloaded = true;
		//		eAction.countRestricted = restrictedAccessCount;
		Globals.session().saveOrUpdate(eAction);
	}

	@Override
	public void setStatus(String newStatus)
	{
		if (eAction != null)
		{
			eAction.status = newStatus;
			Globals.session().saveOrUpdate(eAction);
		}
		super.setStatus(newStatus);
	}
}
