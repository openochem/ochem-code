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

import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.dao.Repository;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.Molecule;
import qspr.frontend.WebList;
import qspr.frontend.WebModel;
import qspr.workflow.utils.QSPRConstants
;
import qspr.modelling.ModelStatistics;
import qspr.modelling.PointStatistics;
import qspr.modelling.SetStatistics;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.applier.ModelApplierTaskProcessor;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.business.PredictionNeighbor;
import com.eadmet.business.PredictionSpaceKnnCalculator;

@Controller
public class ModelNeighboursController extends BrowserWrapper
{

	public ModelNeighboursController()
	{
		sessionRequired = true;
	}
	
	public ModelAndView show(HttpServletRequest request, HttpServletResponse response)
	{
		return new WebModel().setTemplate("model/neighbours").getModelAndView();
	}
	
	public ModelAndView list(HttpServletRequest request, HttpServletResponse response) throws IOException
	{
		List<ExperimentalProperty> eps = new ArrayList<ExperimentalProperty>();
		
		Long model_id = getLongParam("model_id");
		Long mm_id = getLongParam("mm_id");
		
		Model model = (Model) Globals.session().get(Model.class, model_id);
		ModelMapping selectedMapping = model.getMappingById(mm_id);
		ModelStatistics ms = (ModelStatistics) selectedMapping.statisticsRecalculated.getObject();
		SetStatistics trainingSet = ms.sets.get(0);
		List<PredictionNeighbor> knn = new ArrayList<PredictionNeighbor>();
		
		String type = "substructure";
		if (assertParam("type"))
			type = getParam("type");
		
		if (assertParam("ep_id"))
		{
			Long ep_id = getLongParam("ep_id");
			PointStatistics originalPoint = ms.getPointById(ep_id);
			if (!type.equals("structural-similarity"))
			{
				knn = PredictionSpaceKnnCalculator.getNeighbours(originalPoint, trainingSet, type);
			} else
			{
				ExperimentalProperty ep = (ExperimentalProperty)Globals.session().get(ExperimentalProperty.class, originalPoint.id);
				knn = PredictionSpaceKnnCalculator.getStructuralNeighbours(trainingSet, ep.molecule.getData(),ep.property.id);
			}
			
			if (knn.size() > 0)
			{
				int index = 0;
				while (index < knn.size())
					if (knn.get(index).trainingSetPoint.id == originalPoint.id)
						knn.remove(index);
					else
						index++;
			}
			
		} else
		if (assertParam("task_num") && assertParam("row_num"))
		{
			Integer task_num = getIntParam("task_num");
			Integer row_num = getIntParam("row_num");
			
			ModelApplier applier = (ModelApplier) Globals.getSessionAttribute(SessionVariable.MODEL_APPLIER);
			
			if (!type.equals("structural-similarity"))
			{
				ModelApplierTaskProcessor mt = applier.modelTasks.get(task_num);
				DataTable modelResult = mt.wndResult.ports.get(0);
				
				int arrIndex = modelResult.getColumnIndex(QSPRConstants.INDIVIDUAL_PREDICTIONS);
				if (model.modelMappings.size() > 1)
					arrIndex = modelResult.getColumnIndex(QSPRConstants.INDIVIDUAL_PREDICTIONS + selectedMapping._class);
				
				if (arrIndex != -1)
				{
					float[] ensemblePredictions = (float[])modelResult.getRow(row_num).getValue(arrIndex);
					knn = PredictionSpaceKnnCalculator.getNeighbours(ensemblePredictions, trainingSet, type);
				}
			} else
			{
				Long molId = applier.compoundsProvider.basket.entries.get(row_num).ep.molecule.id;
				Molecule m = Repository.molecule.getMolecule(molId);
				knn = PredictionSpaceKnnCalculator.getStructuralNeighbours(trainingSet, m.getData(),null);
			}
		}
		
		for (PredictionNeighbor neighbour : knn.subList(0, Math.min(10, knn.size()))) 
		{
			ExperimentalProperty ep = (ExperimentalProperty) Globals.session().get(ExperimentalProperty.class, neighbour.trainingSetPoint.id);
			neighbour.trainingSetPoint.transferDataToEP(ep, model, selectedMapping, trainingSet);

			DecimalFormat formatter = (DecimalFormat)NumberFormat.getInstance();
			formatter.applyPattern("###,##0.00");
			ep.correl = formatter.format(neighbour.distance);
			
			eps.add(ep);
		}
		
		
		WebList wl = new WebList();
		wl.pageNum = 0;
		wl.pageSize = 10;
		wl.loadFromList(eps);
		
		return new WebModel(wl).getModelAndView();
	}
}