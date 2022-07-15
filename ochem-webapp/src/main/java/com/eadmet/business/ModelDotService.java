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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.entities.BasketEntry;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.workflow.utils.QSPRConstants
;
import qspr.modelling.ModelStatistics;
import qspr.modelling.PointStatistics;
import qspr.modelling.SetStatistics;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.NumericalValueStandardizer;

public class ModelDotService 
{
	public static List<Long> getErrorRecordsIds(ModelMapping modelMapping, boolean recalculated, String errorFilter)
	{
		List<Long> idList = new ArrayList<Long>();

		ModelStatistics modelStat = recalculated ? 
				(ModelStatistics)modelMapping.statisticsRecalculated.getObject() : 
					(ModelStatistics)modelMapping.statisticsOriginal.getObject();
				for (int set = 0; set < modelStat.sets.size(); set++)
					for (int i = 0; i < modelStat.sets.get(set).points.size(); i++)
					{
						PointStatistics pointStat = modelStat.sets.get(set).points.get(i);
						if (pointStat.error != null)
							if (errorFilter == null || pointStat.error.toLowerCase().contains(errorFilter.toLowerCase()))
								idList.add(pointStat.id);
					}

				return idList;
	}

	//An ugly and non-reusable method, since it returns string and not proper name/value pairs. Consider rewriting. Alson consider including pagination filter.
	public List<String> getDescriptorList(Model model, Long modelMappingId, Long epId)
	{
		List<String> descriptorList = new ArrayList<String>();
		if (model != null)
		{
			if (model.calcDescriptors != null)
			{
				ModelMapping mm = model.getMappingById(modelMappingId);
				ModelStatistics stats = ModelStatistics.get(mm);
				DataTable calculatedDescriptors = model.getCalculatedDescriptors();
				calculatedDescriptors.reset();


				if (calculatedDescriptors.getRow(0).attachments.containsKey(QSPRConstants.RECORD_ID_ATTACHMENT))
				{
					// Find the descriptor row by the record ID
					// In future, may be optimize it by hashing record IDs
					// Midnighter
					do 
						try
					{
							calculatedDescriptors.forceNextRow();
					}
					catch (Exception e)
					{
						throw new UserFriendlyException("The descriptors for this point are unavailable");
					}
					while 
						(!epId.equals(calculatedDescriptors.getCurrentRow().getAttachment(QSPRConstants.RECORD_ID_ATTACHMENT)));
				}
				else
				{
					// Backwards compatibility for old models, which do not store record IDs
					// This is the old non-stable logic, does not work properly for multi-learning models

					PointStatistics ps = stats.getPointById(epId);
					calculatedDescriptors.currentRow = ps.numInSet;
					if (ps.parent.setId.equals(QSPRConstants.VALIDATION))
						calculatedDescriptors.currentRow += stats.sets.get(0).points.size();
				}

				for (String column : calculatedDescriptors.getColumns())
				{
					Serializable desc = calculatedDescriptors.getValue(column);

					if (!(desc instanceof Double))
						continue;
					descriptorList.add(column+" = "+NumericalValueStandardizer.getSignificantDigitsStr((Double) desc, 3));
				}
			}				
		}
		return descriptorList;
	}

	@SuppressWarnings("unchecked")
	public List<ExperimentalProperty> getErrorList(Model model, PaginationFilter pager, boolean recalculated, String errorTitle)
	{
		Map<Long, String> errorMap = new HashMap<Long,String>();
		List<ExperimentalProperty> globalList = new ArrayList<ExperimentalProperty>();

		pager.totalSize = 0L;

		List<ModelMapping> modelList = model.modelMappings;
		for (ModelMapping modelMapping : modelList) 
		{
			int totalErrors = 0;
			ModelStatistics modelStat = recalculated ?
					(ModelStatistics)modelMapping.statisticsRecalculated.getObject() : 
						(ModelStatistics)modelMapping.statisticsOriginal.getObject();

					List<Long> idList = new ArrayList<Long>();

					for (int set = 0; set < modelStat.sets.size(); set++)
						for (int i=0; i < modelStat.sets.get(set).points.size(); i++)
						{
							PointStatistics pointStat = modelStat.sets.get(set).points.get(i);
							if (pointStat.error != null)
								if (errorTitle == null || pointStat.error.toLowerCase().contains(errorTitle.toLowerCase()))
								{
									totalErrors++;
									errorMap.put(pointStat.id, pointStat.error);
									if (pager.pageSize == 0 || (totalErrors >= (pager.pageNum - 1)*pager.pageSize + 1 && totalErrors <= (pager.pageNum*pager.pageSize)))
									{
										idList.add(pointStat.id);
									}
								}
						}

					List<ExperimentalProperty> recordList  = null;

					if (idList.size() > 0)
					{
						recordList = Globals.session().createCriteria(ExperimentalProperty.class)
								.add(Restrictions.isNull("deleted"))
								.add(Restrictions.in("id", idList)).list();

						for (ExperimentalProperty record : recordList) 
						{
							record.error = errorMap.get(record.id);
							globalList.add(record);
						}
					}

					pager.totalSize += totalErrors;
		}

		return globalList;

	}

	@SuppressWarnings("unchecked")
	public List<ExperimentalProperty> getModelDotList(Model model, Long statNum, Long epId, Long modelMappingId)
	{
		ExperimentalProperty ep = null;
		List<ExperimentalProperty> eps = new ArrayList<ExperimentalProperty>();	
		if (model != null && model.attachment != null)
		{
			ModelStatistics ms; //= (ModelStatistics) model.getAttachment().get("statisticsOriginal");
			//long setNum = getLongParam("setnum");
			//long pointNum = getLongParam("pointnum");

			//TODO: Make sure original descriptors are stored somewhere

			ModelMapping selectedMapping = model.getMappingById(modelMappingId);
			if (statNum == 0)
				ms = (ModelStatistics) selectedMapping.statisticsOriginal.getObject();
			else
				ms = (ModelStatistics) selectedMapping.statisticsRecalculated.getObject();
			//ms.actualizeStatistics(selectedMapping);
			//SetStatistics ss = ms.sets.get((int) setNum);
			PointStatistics originalPoint = ms.getPointById(epId); //ss.points.get((int)pointNum);
			long moleculeId = originalPoint.getMoleculeId();
			long recordId = originalPoint.id;

			List<ModelMapping> modelList = model.modelMappings;

			for (ModelMapping modelMapping : modelList) 
			{
				if (statNum == 0 || modelMapping.statisticsRecalculated == null)
					ms = (ModelStatistics) modelMapping.statisticsOriginal.getObject();
				else
					ms = (ModelStatistics) modelMapping.statisticsRecalculated.getObject();

				for (SetStatistics ss : ms.sets)
				{
					//SetStatistics ss = originalPoint.parent;

					for (PointStatistics ps : ss.points) 
					{
						if ((ps.getMoleculeId() != moleculeId) || (moleculeId == 1 && ps.id != recordId))
							continue;

						ep = (ExperimentalProperty)Globals.session().get(ExperimentalProperty.class, ps.id);

						if (ep != null && ep.deleted == null)
						{
							ps.transferDataToEP(ep, model, modelMapping, ss);//??
							ep.descriptorList = new DataTable();
							//add exclude/include information
							List <BasketEntry> beList = 
									Globals.session().createCriteria(BasketEntry.class)
									.add(Restrictions.eq("basket", model.trainingSet))
									.createCriteria("ep").add(Restrictions.eq("id", ep.id)).list();

							if (beList.size() > 0)
								ep.exclude = model.microattachment.getObject().excludedBasketEntries.contains(beList.get(0).id) || (beList.get(0).exclude);

							//TODO: Put check like "neighbours exist" or something

							if (ps == originalPoint)
								eps.add(0, ep); // The clicked point goes first in the list
							else
								eps.add(ep);

						}
					}
				}
			}
		}
		return eps;
	}
}

