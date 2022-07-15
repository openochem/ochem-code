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

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.hibernate.Hibernate;
import org.hibernate.criterion.Projections;
import org.hibernate.criterion.Restrictions;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Attachment;
import qspr.entities.Attachment.AttachmentType;
import qspr.entities.AttachmentSource;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.DataHandlingOptions;
import qspr.entities.ExperimentalProperty;
import qspr.entities.Model;
import qspr.entities.ModelAttachment;
import qspr.entities.ModelMapping;
import qspr.entities.ModelMicroAttachment;
import qspr.entities.Molecule;
import qspr.entities.Property;
import qspr.entities.PropertyOption;
import qspr.entities.Unit;
import qspr.metaserver.configurations.ASNNConfiguration;
import qspr.metaserver.configurations.SelectionConfiguration;
import qspr.workflow.utils.QSPRConstants
;
import qspr.modelling.ModelUploadProcessor;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.util.BasicRecordMapper;
import qspr.util.MoleculePeer;
import qspr.util.StatusTracker;
import qspr.util.UploadContext;
import qspr.workflow.datatypes.CompactDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.exceptions.UserFriendlyException;
import com.eadmet.parsers.SimpleParser;

public class ModelUploadService 
{
	public StatusTracker status = new StatusTracker();

	private Model prepareModelStub(Long trainingSetId, List<Long> validationSetIds, Map<String, Long> unitIds)
	{
		Model model = new Model();
		model.session = Globals.userSession();

		if (trainingSetId != null)
			model.trainingSet = Basket.getBasket(Globals.userSession(), trainingSetId);

		model.cleanValidationSets();
		if (validationSetIds != null && !validationSetIds.isEmpty())
			for (Long valSetID : validationSetIds)
				model.addValidationSet(Basket.getBasket(Globals.userSession(), valSetID));

		BasicRecordMapper brm = new BasicRecordMapper(model.trainingSet);
		List<Property> propList = model.trainingSet.getProperty();
		int count = 0;
		model.modelMappings.clear();
		for (Property property : propList.subList(0, 1))
		{
			Hibernate.initialize(property.options);
			ModelMapping modelMapping = new ModelMapping();
			modelMapping.property = property;

			if (unitIds.containsKey("unit"+count))
				modelMapping.unit = (Unit) Globals.session().get(Unit.class, unitIds.get("unit"+count));
			else
				modelMapping.unit = property.defaultUnit;
			modelMapping._class = brm.getClass(modelMapping.property);
			model.addModelMapping(modelMapping);
			count++;
		}

		model.attachment = new Attachment<ModelAttachment>(new ModelAttachment(), AttachmentType.MARSHALABLE, AttachmentSource.Model);
		model.microattachment = new Attachment<ModelMicroAttachment>(new ModelMicroAttachment(), AttachmentType.MARSHALABLE, AttachmentSource.Model);

		model.attachment.getObject().configuration =  new CDSConfiguration();

		model.attachment.updateObject();
		model.template = Repository.modelTemplate.getByName(QSPRConstants.UPLOADED_MODEL);
		model.name = QSPRConstants.UPLOADED_MODEL;
		return model;
	}

	@SuppressWarnings("unchecked")
	private void fillMapsFromBaskets(List<Long> involvedEpIds, Map<String, List<Long>> externalIdMap, Map<Integer, List<Long>> mappingIdMap, List<Long> baskets)
	{
		List<Object[]> list = 
				Globals.session().createCriteria(ExperimentalProperty.class)
				.createAlias("basketEntries", "be")
				.createAlias("be.basket", "b")
				.createAlias("molecule", "m")
				.createAlias("m.mapping2", "m2")
				.add(Restrictions.in("b.id", baskets))
				.setProjection(
						Projections.projectionList()
						.add(Projections.property("id"))
						.add(Projections.property("m2.id"))
						.add(Projections.property("externalId")))
						.list();

		for (Object[] res : list) 
		{
			Long epId = (Long)res[0];
			Integer mpId = (Integer)res[1];
			String extId = (String)res[2];

			involvedEpIds.add(epId);

			if (extId != null)
			{
				if (externalIdMap.get(extId) == null)
					externalIdMap.put(extId, new ArrayList<Long>());
				externalIdMap.get(extId).add(epId);
			}

			if (mpId != null)
			{
				if (mappingIdMap.get(mpId) == null)
					mappingIdMap.put(mpId, new ArrayList<Long>());
				mappingIdMap.get(mpId).add(epId);
			}
		}
	}

	private WorkflowNodeData getModelWnd(ModelUploadResult result, File f) throws Exception
	{
		Model model = result.model;
		String propertyName = model.modelMappings.get(0).property.getName();
		DataTable dtTmpPredictions = new DataTable(true);
		dtTmpPredictions.addColumn(QSPRConstants.PREDICTION_RESULT_COLUMN);


		List<Long> baskets = new ArrayList<Long>();
		baskets.add(model.trainingSet.id);
		for (Basket b : model.getValidationSets()) 
			baskets.add(b.id);

		List<Long> involvedEpIds = new ArrayList<Long>();
		Map<String, List<Long>> externalIdMap = new HashMap<String, List<Long>>();
		Map<Integer, List<Long>> mappingIdMap = new HashMap<Integer, List<Long>>();

		fillMapsFromBaskets(involvedEpIds, externalIdMap, mappingIdMap, baskets);

		Map<String, Long> propertyOptionsMap = new HashMap<String, Long>();
		List<PropertyOption> lpo = model.modelMappings.get(0).property.options;
		for (PropertyOption propertyOption : lpo)
			propertyOptionsMap.put(propertyOption.name, propertyOption.id);


		SimpleParser parser = SimpleParser.getParser(f.getAbsolutePath()).setSource(f);
		parser.setCurrentSheet(0);
		List<String> headers = parser.sheetColumns.get(0);

		int extIdIndex = -1; 
		int mpIdIndex = -1;
		int predIndex = -1;

		for (int i=0; i<headers.size(); i++) 
			if (headers.get(i).equalsIgnoreCase("external_id"))
				extIdIndex = i;
			else if (headers.get(i).equalsIgnoreCase("molecule") || headers.get(i).equalsIgnoreCase("molecule_id") || headers.get(i).equalsIgnoreCase("moleculeid"))
				mpIdIndex = i;
			else if (headers.get(i).equalsIgnoreCase("prediction") || headers.get(i).matches(propertyName+" .*"))
				predIndex = i;

		if (mpIdIndex == -1 && extIdIndex == -1)
		{
			result.errors.add("Neither MOLECULE nor EXTERNAL_ID columns found in uploaded file");
			return null;
		}

		if (predIndex == -1)
		{
			result.errors.add("PREDICTION column not found in uploaded file");
			return null;
		}

		int rowNum = 1;

		Map<Long, Integer> epToRow = new HashMap<Long, Integer>();

		Set<Long> foundEP = new HashSet<Long>();

		for (List<String> row : parser) 
		{
			String prediction = row.get(predIndex);
			List<Long> expPropertyId = null;

			if (extIdIndex != -1)
			{
				String extId = row.get(extIdIndex);
				if (externalIdMap.get(extId) == null){
					result.errors.add("Row "+rowNum+": Could not find existing experimental record in training or validation sets for external id " + extId);
					continue;
				}

				expPropertyId = externalIdMap.get(extId);
			} else
				if (mpIdIndex != -1)
				{
					String mol = row.get(mpIdIndex);
					try
					{
						Molecule m = MoleculePeer.fetchFromString(mol,  new UploadContext());

						if (mappingIdMap.get(m.mapping2.id) == null){
							// nor experimental record, skipping ...
//							result.errors.add("Row "+rowNum+": Could not find existing experimental record in training or validation sets for molecule "+mol);
							continue;
						}

						expPropertyId = mappingIdMap.get(m.mapping2.id);
					} catch (Exception e)
					{
						result.errors.add("Row "+rowNum+": "+e.getMessage());
						continue;
					}
				}

			if (expPropertyId == null){
				result.errors.add("Row "+rowNum+": Could not match uploaded predicted value to the existing experimental record");
				continue;
			}

			rowNum++;

			for(Long id : expPropertyId){
				if(foundEP.contains(id))continue; // already added; it is a duplicate
				foundEP.add(id);

				dtTmpPredictions.addRow();
				dtTmpPredictions.getCurrentRow().addAttachment(QSPRConstants.RECORD_ID_ATTACHMENT, id);

				if (model.modelMappings.get(0).property.isNumeric())
					dtTmpPredictions.setValue(Double.valueOf(prediction));
				else
				{
					Long optionId = propertyOptionsMap.get(prediction);
					if (optionId == null)
					{
						result.errors.add("Row "+rowNum+": Unknown option "+prediction+" for qualitative property "+model.modelMappings.get(0).property.getName());
						expPropertyId = null;
					}
					else
						dtTmpPredictions.setValue(Double.valueOf(model.getMappedOption(optionId)));
				}

				epToRow.put(id, dtTmpPredictions.currentRow);
			}
		}

		while(involvedEpIds.removeAll(foundEP)); // can be required to remove several elements two times

	/*	// we allow to have not uploaded predictions for some records
		if (involvedEpIds.size() > 0)
		{
			StringBuilder error = new StringBuilder();
			error.append("The following experimental records from the training and validation sets did not have predicted values in the uploaded file: ");
			for (Long epid : involvedEpIds)
				error.append("R"+epid+", ");
			result.errors.add(error.substring(0, error.length() - 2));
		}
	*/	
		if (result.errors.size() > 0)
			return null;

		DataTable dtPredictions = new DataTable(true);
		dtPredictions.addColumn(QSPRConstants.PREDICTION_RESULT_COLUMN);

		for (BasketEntry be : model.trainingSet.entries)
			if(epToRow.get(be.ep.id) != null)dtPredictions.addRow(dtTmpPredictions.getRow(epToRow.get(be.ep.id)));
			else
				dtPredictions.addRow((new CompactDataRow()).setError("No value was uploaded for this record.").addAttachment(QSPRConstants.RECORD_ID_ATTACHMENT, be.ep.id));

		for (Basket vs : model.getValidationSets())
			for (BasketEntry be : vs.entries)
				if(epToRow.get(be.ep.id) != null)dtPredictions.addRow(dtTmpPredictions.getRow(epToRow.get(be.ep.id)));
				else
					dtPredictions.addRow((new CompactDataRow()).setError("No value was uploaded for this record.").addAttachment(QSPRConstants.RECORD_ID_ATTACHMENT, be.ep.id));

		WorkflowNodeData uploadedWnd = new WorkflowNodeData();
		uploadedWnd.addPort(dtPredictions);
		uploadedWnd.addPort(new DataTable (new ASNNConfiguration()));
		uploadedWnd.addPort(dtPredictions.getCopy().setId("descriptors"));
		uploadedWnd.addPort(new DataTable(new SelectionConfiguration()).setId("selection-configuration"));

		return uploadedWnd;
	}

	public ModelUploadResult uploadModel(ModelUploadQuery query) throws Exception
	{
		if(query.file == null)throw new UserFriendlyException("The file with predicted data has not been provided or found. Please, check the file and try again.");

		ModelUploadResult result = new ModelUploadResult();

		status.set("Preparing model stub");

		Model model = result.model = prepareModelStub(query.trainingSetId, query.validationSetIds, query.unitMap);

		ModelUploadProcessor processor = new ModelUploadProcessor();
		processor.statusTracker.addListener(status);

		processor.model = model;
		processor.prepare();

		DataHandlingOptions dho = model.attachment.getObject().datahandling;
		dho.approximateequals = QSPRConstants.USE;
		dho.greaterless = QSPRConstants.USE;
		dho.intervals = QSPRConstants.USE;
		model.attachment.updateObject();

		WorkflowNodeData uploadedWnd = getModelWnd(result, query.file);

		if (uploadedWnd == null) //Some errors
			return result;

		processor.setUploadedData(uploadedWnd);
		processor.start();
		processor.join(0);

		status.set("Saving uploaded model");

		Globals.session().delete(processor.pTask);
		model.attachment.getObject().configuration = null;
		model.attachment.getObject().datahandling = null;
		model.attachment.getObject().protocol = null;
		model.attachment.getObject().standartization = null;
		model.attachment.updateObject();
		model.description = QSPRConstants.UPLOADED_MODEL;
		model.userDescription = query.description;
		model.configurationXml = "";
		model.configurationHash = "";
		model.taskId = null;
		model.name = query.description + " file="+query.file.getName(); 
		Globals.session().saveOrUpdate(model);

		return result;
	}


}
