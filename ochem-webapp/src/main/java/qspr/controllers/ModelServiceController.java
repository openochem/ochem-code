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

import java.util.HashMap;
import java.util.Map;

import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.stereotype.Controller;
import org.springframework.web.servlet.ModelAndView;
import org.springframework.web.servlet.mvc.multiaction.MultiActionController;

import qspr.dao.Various;
import qspr.entities.ModelIdentity;
import qspr.services.ModelResponse;
import qspr.services.ModelService;
import qspr.services.ModelSummary;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.utils.OCHEMUtils;

/**
 * This is a simple wrapper around SOAP web services
 * that allows to run predictions via a much simpler REST interface
 * and get prediction results in JSON format 
 * 
 * @author midnighter
 *
 */
@Controller
public class ModelServiceController extends MultiActionController
{
	private static transient final Logger logger = LogManager.getLogger(ModelServiceController.class);

	/**
	 * An in-memory cache of the tasks. The key is an MD5 of the requested model ID and the molecules.
	 */
	private static Map<String, ModelResponse> resultsCache = new HashMap<String, ModelResponse>();
	private static Map<String, ModelAndView> taskCache = new HashMap<String, ModelAndView>();


	private ModelAndView getResponsePrediction(String session,long modelId, String sdf, boolean fetchAlsoResults) {
		ModelResponse result = null;

		String inputMD5 = OCHEMUtils.getMD5(sdf) + modelId;

		if (resultsCache.get(inputMD5) == null)
		{
			logger.info("New calculations will start for session: " + session + " mols: " + sdf);

			String sdfs[] = OCHEMUtils.splitSDF(sdf);

			for(int i =0; i<sdfs.length; i++)try {
				sdfs[i] = Various.molecule.convertToKekuleSMILES(Various.molecule.convertToCanonicalSMILES(sdfs[i]));
			}catch(Exception e) {
			}

			result = new ModelService().postModelWithSession(session, modelId, sdfs);
			if(!fetchAlsoResults || result.getStatus().equals(QSPRConstants.ERROR_STATUS)) // do we need to process further?
				return new ModelAndView("json", "object", result);
			resultsCache.put(inputMD5, result);
			result = new ModelService().fetchModel(session, result.getTaskId());
			if(result.getStatus().equals("success"))
				resultsCache.put(inputMD5, result);
		}else
			logger.info("A cached result was found for session: " + session);

		result = resultsCache.get(inputMD5);

		if(result.getTaskId() != 0) {
			result = new ModelService().fetchModel(session, result.getTaskId());
			if(result.getStatus().equals("success") || result.getStatus().equals("failed"))
				resultsCache.put(inputMD5, result);
		}

		logger.info("Providing results for session: " + session + " mols: " + sdf);

		return new ModelAndView("json", "object", result);

	}

	/**
	 * Run predictions with authentification
	 */
	public ModelAndView getPredictionAuth(HttpServletRequest request, HttpServletResponse response) {

		String session = request.getParameter("session");
		if(session == null) return noSessionProvided();
		logger.info("starting a new request with session: " + session);
		return getResponsePrediction(session, getModelID(request), request.getParameter("mol"),true);
	}

	/**
	 * Run predictions in a simplified way, without bothering about task IDs 
	 */
	public ModelAndView getPrediction(HttpServletRequest request, HttpServletResponse response) {
		return getResponsePrediction(QSPRConstants.ANONYMOUS, getModelID(request), request.getParameter("mol"),true);
	}

	public ModelAndView getModelDescriptionAuth(HttpServletRequest request, HttpServletResponse response) {
		long modelId = getModelID(request);
		String session = request.getParameter("session");
		if(session == null) return noSessionProvided();
		ModelSummary[] summary = new ModelService().getModelSummary(session, modelId);
		return new ModelAndView("json", "object", summary);
	}

	ModelAndView noSessionProvided() {
		ModelResponse result = new ModelResponse();
		result.setStatus(QSPRConstants.ERROR_STATUS);
		result.setDetailedStatus(QSPRConstants.NOSESSION + " provide session id to run calculations using this service");
		return new ModelAndView("json", "object", result);
	}


	private long getModelID(HttpServletRequest request) {
		String modelId = request.getParameter("modelId");
		logger.info("Posting a model " + modelId);
		Long publicId = null;
		if (modelId.matches("[0-9]+"))
			publicId = Long.valueOf(modelId);
		else
		{
			ModelIdentity mi = ModelIdentity.getByGUID(modelId);
			if (mi != null)
				publicId = mi.model.publicId;
		}

		if (publicId == null)
			throw new UserFriendlyException("Unknown model identity: " + modelId);

		return publicId;
	}

	// will not work without session


	/**
	 * Another way to run predictions using task IDs (kept for compatibility)
	 * Required for model services
	 */
	public ModelAndView postModel(HttpServletRequest request, HttpServletResponse response)
	{
		long modelId = getModelID(request);
		String sdf = request.getParameter("mol");
		String inputMD5 = OCHEMUtils.getMD5(sdf) + modelId;
		if (taskCache.get(inputMD5) == null)
			taskCache.put(inputMD5, getResponsePrediction(null, modelId, sdf,false));

		return taskCache.get(inputMD5);
	}


	public ModelAndView fetchModel(HttpServletRequest request, HttpServletResponse response)
	{
		String session = request.getParameter("session");
		if(session == null) return noSessionProvided();

		logger.info("Fetching task " + request.getParameter("taskId"));
		ModelResponse result = new ModelService().fetchModel(session, Long.valueOf(request.getParameter("taskId")));
		ModelAndView mav = new ModelAndView("json", "object", result);
		return mav;
	}

}
