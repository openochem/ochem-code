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

import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;

import qspr.Globals;
import qspr.SessionVariable;
import qspr.dao.Repository;
import qspr.frontend.WebModel;
import qspr.modelling.CompoundsProvider;
import qspr.modelling.applier.ModelApplier;



// Midnighter on Jun 26, 2012
/**
 * A simplified dialog that allows to run the best of our models
 * This dialog is focused on the predicted endpoints rather than on models
 * It is going simplify user's life a lot
 * 
 * @author midnighter
 */
@Controller
public class PredictorController extends ControllerWrapper
{
	public PredictorController()
	{
		this.sessionRequired = true;
	}
	
	public ModelAndView show(HttpServletRequest request, HttpServletResponse response)
	{
		return new WebModel().setList(Repository.model.getFeaturedModels()).setTemplate("model/predictor").getModelAndView();
	}
	
	public ModelAndView start(HttpServletRequest request, HttpServletResponse response)
	{
		ModelApplier applier = new ModelApplier();
		Globals.setSessionAttribute(SessionVariable.MODEL_APPLIER, applier);
		String[] modelIDsStr = request.getParameterValues("model");
		
		for (String modelIDStr : modelIDsStr) {
			applier.addModel(Repository.model.getByPublicId(Long.valueOf(modelIDStr)));
		}
		
		if (assertParam("disable-cache"))
			applier.useCache = false;
		
		CompoundsProvider cmpProvider = new CompoundsProvider();
		applier.compoundsProvider = cmpProvider.parseUI();
		
		return new WebModel().setTemplate("model/apply").getModelAndView();
	}
}
