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

package qspr.metaserver.cs.util;

import java.io.IOException;

import ochem.eadmet.wsapi.DataReferenceRequest;
import ochem.eadmet.wsapi.ModelResponse;
import ochem.eadmet.wsapi.ModelServicePortType;
import ochem.eadmet.wsapi.Prediction;
import ochem.eadmet.wsapi.PropertyPrediction;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.metaserver.configurations.ModelAbstractConfiguration.PredictionScenario;
import qspr.metaserver.configurations.DescriptorType;
import qspr.metaserver.configurations.DescriptorsApplyModelConfiguration;
import qspr.metaserver.cs.WorkflowNodeServer;
import qspr.metaserver.protocol.Task;
import qspr.metaserver.transport.DataReference;
import qspr.metaserver.transport.DataReferenceFactory;
import qspr.metaserver.transport.DataReferencer;
import qspr.metaserver.util.ServiceLocator;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.exceptions.UserFriendlyException;

public class WebServiceCalculationTask extends DescriptorsCacheTask
{

	protected static final Logger logger = LogManager.getLogger(WebServiceCalculationTask.class);

	private Long parentTaskId;
	private int countFailure =0;
	Boolean useCache;
	private String sessionId;

	static private WorkflowNodeData savedWndInput;
	static private String savedNoSQLReference;

	static final int MAXIMUM_FAILURES = 10;

	ModelServicePortType service;
	public ModelResponse response;

	public WebServiceCalculationTask(DescriptorType deskType,
			WorkflowNodeData wndinput,
			WorkflowNodeServer theParent,Boolean forceUpdateCache) throws Exception {
		super(deskType,wndinput,theParent,forceUpdateCache);
	}


	@Override
	public void post() throws IOException
	{
		if(wndOutput != null)return; // all data are found in cache; no posting is required

		logger.info("Starting posting " + taskName+ " " +configuration);

		DescriptorsApplyModelConfiguration configuration = (DescriptorsApplyModelConfiguration) this.configuration;

		sessionId = (String)wndInput.ports.get(0).getColumnAttachment(QSPRConstants.SDF_COLUMN, "session-guid");

		if(wndInput != savedWndInput){ // in case we have modified data because of cache
			DataReferencer referencer = DataReferenceFactory.createReferencer();
			DataReference ref = referencer.saveReference(wndInput, QSPRConstants.WEBSERVICE_DATABASE);
			logger.info("Data size in " + referencer.getDataSize(ref));
			savedNoSQLReference = ref.getReference();
			savedWndInput = wndInput;
		}

		DataReferenceRequest  predictionRequest = new DataReferenceRequest();
		predictionRequest.setModelId(configuration.modelId);
		predictionRequest.setSessionGUID(sessionId);
		predictionRequest.setDataReference(savedNoSQLReference);
		predictionRequest.setDatasize(wndInput.ports.get(0).getRowsSize());
		if (configuration.scenario != null)
			predictionRequest.setPredictionScenario(configuration.scenario.toString());
		else
			predictionRequest.setPredictionScenario(PredictionScenario.PREDICTION_ONLY.toString());

		for(int i=0;i<=MAXIMUM_FAILURES;i++)
			try
		{
				Thread.sleep(i*100*1000);
				if(service == null) service = ServiceLocator.getModelService();
				response = service.postDataReferenceRequest(predictionRequest);
				if (response.getStatus().equals(QSPRConstants.SUCCESS))
				{
					taskId = response.getMetaserverTaskId().intValue(); // this is in the Metaserver Database
					parentTaskId = response.getTaskId(); // this is in the Pending Tasks
					Task t = parent.currentTask.get();
					int parentId = 0;

					if(t != null)
						parentId = t.id > 0?t.id: t.parentTaskId != null && t.parentTaskId > 0? t.parentTaskId :0;

						if (parentId > 0){
							client.setTaskParent(taskId, parentId);
							client.setTaskPriority(taskId, (int) Math.round(t.priority)); // correcting priority
						}
						return;
				}
				else
					throw new UserFriendlyException(response.getDetailedStatus());
		}
		catch (UserFriendlyException e)
		{
			throw e;
		}
		catch (Exception e)
		{
			logger.info("Could not post the " + taskName+ " "+ configuration+" trial="+i);

			e.printStackTrace();
			if(++i>MAXIMUM_FAILURES)
				throw new UserFriendlyException("Could not post the " + taskName+ " "+ configuration);
		}
	}

	@Override
	public void cleanup()
	{
		service = null;
		savedWndInput = null;
		savedNoSQLReference = null;
		System.gc();
	}

	@Override
	public boolean isReady() throws Exception
	{
		if(wndOutput != null || Task.ERROR.equals(status))return true; // work is done!

		try
		{
			DescriptorsApplyModelConfiguration configuration = (DescriptorsApplyModelConfiguration) this.configuration;
			response = service.getPredictions(sessionId, parentTaskId, configuration.units);
			status = response.getDetailedStatus();

			if (response.getStatus() != null && response.getStatus().equals("pending"))
			{
				countFailure=0; 
				return false;
			}

			if (response.getStatus() != null && response.getStatus().equals(Task.ERROR)){
				status = Task.ERROR;
				return true;
			}

			DataTable dtPredictions = new DataTable(true);
			for (Prediction prediction : response.getPredictions())
			{
				dtPredictions.addRow();
				if (prediction.getError() != null)
					dtPredictions.getCurrentRow().setError(prediction.getError());
				else
					for (PropertyPrediction pp : prediction.getPredictions()){ // storing both predictions and their accuracy
						dtPredictions.setValue(pp.getProperty() + " (model "+((DescriptorsApplyModelConfiguration) configuration).modelId+")", pp.getValue());
						dtPredictions.setValue(pp.getProperty() + "-accuracy (model "+((DescriptorsApplyModelConfiguration) configuration).modelId+")", pp.getAccuracy());
					}
			}
			wndOutput = new WorkflowNodeData(dtPredictions);
			service.deleteTask(sessionId, parentTaskId);

			logger.info("Finished " + taskName+ " " +configuration + " " + dtPredictions.getRowsSize() + 
					(dtPredictions.getRowsNoErrorsSize() != dtPredictions.getRowsSize()? " errors: " + 
							(dtPredictions.getRowsSize() - dtPredictions.getRowsNoErrorsSize()) : ""));

			onReady(); // add previously cached or failed predictions

			return true;
		}	
		catch (Exception e) {
			logger.error(e);
			if(++countFailure>MAXIMUM_FAILURES){
				logger.info("Finished " + taskName+ " " +configuration + " - calculation failed");
				throw new UserFriendlyException("Failed to retrieve predictions for " + taskName + " " +
						configuration + " because of error: "+e.getLocalizedMessage());
			}
			return false;
		}
	}
}