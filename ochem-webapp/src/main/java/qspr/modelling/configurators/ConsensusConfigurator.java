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

package qspr.modelling.configurators;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.servlet.http.HttpServletRequest;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import qspr.Globals;
import qspr.business.ModelPeer;
import qspr.dao.ChemInfEngine;
import qspr.dao.Repository;
import qspr.entities.Attachment;
import qspr.entities.AttachmentSource;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.entities.ModelMapping;
import qspr.entities.Unit;
import qspr.exception.UnitConversionException;
import qspr.frontend.WebModel;
import qspr.interfaces.ProvidedConditions;
import qspr.metaserver.configurations.ConsensusModelConfiguration;
import qspr.metaserver.configurations.ConsensusModelConfiguration.ConsensusType;
import qspr.metaserver.configurations.ConsensusModelConfiguration.IndividualModel;
import qspr.workflow.utils.QSPRConstants;
import qspr.metaserver.util.ShortCondition;
import qspr.metaserver.util.aggregator.AveragingAggregator;
import qspr.modelling.ApplicabilityDomain;
import qspr.modelling.ModelStatistics;
import qspr.modelling.PointStatistics;
import qspr.modelling.SetStatistics;
import qspr.util.unitconversion.UnitConversion;
import qspr.workflow.datatypes.DataTable;

import com.eadmet.exceptions.UserFriendlyException;

public class ConsensusConfigurator extends BasicModelConfigurator
{
	private static transient final Logger logger = LogManager.getLogger(ConsensusConfigurator.class);

	public ConsensusConfigurator()
	{
		firstConfigurationStep = "chooseModels";
	}

	@Override
	public void setModel(Model model)
	{
		super.setModel(model);
		model.attachment.getObject().configuration = new ConsensusModelConfiguration();
	}

	public WebModel chooseModels()
	{
		WebModel wm = new WebModel(model);
		if (getDefaultTemplate() != null)
			wm.addObject(getDefaultTemplate());
		return wm.setTemplate("modeller/configurators/choose-consensus-models");
	}

	public void setStatus(String s)
	{
		logger.info(s);
	}

	public void chooseModelsSubmit(HttpServletRequest request) throws Exception
	{
		String[] ids = request.getParameterValues("model-id");
		ConsensusModelConfiguration consConf = (ConsensusModelConfiguration) model.attachment.getObject().configuration;

		consConf.allowErrors = request.getParameter("allow-errors") != null;

		consConf.type = ConsensusType.valueOf(request.getParameter("type"));
		model.attachment.getObject().protocol.validationConfiguration = null;

		for (String id : ids)
			if (!"".equals(id))
			{
				Model model = (Model) Globals.session().get(Model.class, Long.valueOf(id));
				consConf.addModel(model.publicId);
			}

		if (consConf.individualModels.size() < 2)
			throw new UserFriendlyException("A consensus model must include at least two different individual models.");

		createConsensus();
		currentPage = "save";
	}

	/**
	 * Creates the consensus model on the fly, without recalculation of the sub-models
	 * Naturally, this requires compatibility of the sub-models training sets
	 * @throws Exception
	 */

	private void createConsensus() throws Exception 
	{
		ConsensusModelConfiguration consConf = (ConsensusModelConfiguration) model.attachment.getObject().configuration;
		model.name = "Consensus " + ModelPeer.getModelName(model, true);

		Globals.commitAlternateTransaction(); // we do not need them ...

		Repository.model.getByPublicId(consConf.individualModels.get(0).id).restoreModifiedBasketEntries(true); // we always restore the basket

		mergeSets(model);
		logger.info("Downloading individual models...");

		List<Model> childModels = new ArrayList<Model>();
		for (IndividualModel iModel : consConf.individualModels)
			childModels.add(Repository.model.getByPublicId(iModel.id));

		logger.info("Downloaded individual models...");
		
		ApplicabilityDomain[][] ads = new ApplicabilityDomain[childModels.size()][model.modelMappings.size()];
		
		switch(consConf.type) {
		case AVERAGE:
		case RMSE_WEIGHTED:
		case OPTIMAL:
			break;
		case BEST_MODEL:
		case WEIGHTED_AVERAGE:
			logger.info("Checking applicability domains");

			if (!ApplicabilityDomain.hasDM(childModels.get(0)))
				throw new UserFriendlyException("Model " + childModels.get(0) + " does not support prediction accuracy estimates and cannot be a part of an accuracy-weighted consensus. So far, only bagging models or neural network models support prediction accuracy assessment.");

			String dmName = ApplicabilityDomain.getDmName(childModels.get(0));
			for (int i = 0; i < childModels.size(); i++)
				for (int mmNum = 0; mmNum < model.modelMappings.size(); mmNum++)
				{
					if (!ApplicabilityDomain.hasDM(childModels.get(i)))
						throw new UserFriendlyException("Model " + childModels.get(i) + " does not support prediction accuracy estimates and cannot be a part of an accuracy-weighted consensus. So far, only bagging models or neural network models support prediction accuracy assessment.");

					ads[i][mmNum] = new ApplicabilityDomain();
					ads[i][mmNum].setModel(childModels.get(i).modelMappings.get(mmNum), dmName, null);
				}
		}

		for(Model model: childModels)
			if(model.attachment.getObject().configuration instanceof ConsensusModelConfiguration)
				throw new UserFriendlyException("Consensus model should not be used as a part of another consensus. Remove model: " + model);

		model.attachment.getObject().optionsMapping = childModels.get(0).attachment.getObject().optionsMapping;
		model.attachment.getObject().standartization = childModels.get(0).attachment.getObject().standartization;
		
		if(((ProvidedConditions)childModels.get(0).attachment.getObject().configuration).hasConditions()) {
			consConf.conditions = new ArrayList<ShortCondition>();
			for(ShortCondition o: ((ProvidedConditions)childModels.get(0).attachment.getObject().configuration).getConditions()) {
				ShortCondition d = new ShortCondition((ShortCondition)o);
				consConf.conditions.add(d);
			}	
		}

		ChemInfEngine engine = model.attachment.getObject().standartization.synchronise();

		for(Model model: childModels) {
			if(!compareConditions(consConf.conditions,
					((ProvidedConditions)model.attachment.getObject().configuration).getConditions()))
				throw new UserFriendlyException("Can't create a consensus: the model: "+ childModels.get(0) + "\t and model: \t" + model +" have different sets of conditions");
			if(!engine.equals(model.attachment.getObject().standartization.synchronise()))
				throw new UserFriendlyException("Can't create a consensus: the model: "+ childModels.get(0) + "\t and model: \t" + model +" have different sets of standartisation options: "
						+ engine + " and " + model.attachment.getObject().standartization.synchronise());
		}

		if(consConf.type == ConsensusType.RMSE_WEIGHTED) {
			consConf.weights = new ArrayList<Float[]>();
			for	(int i = 0; i < childModels.size(); i++){  // going by individual models for the same property
				Float rmses[] = new Float[model.modelMappings.size()];  // for all properties
				consConf.weights.add(rmses);
			}
		}

		logger.info("Started construction of the consensus model for each property");

		for (int property = 0; property < model.modelMappings.size(); property++) // going by all properties
			selectBestAverage(property, childModels, ads, consConf);
		consConf.compact();

		List<String> columnNames = new ArrayList<String>();

		for(int i=0; i< model.modelMappings.size(); i++){
			String name=model.modelMappings.size()==1?"":""+i;
			columnNames.add(QSPRConstants.PREDICTION_RESULT_COLUMN+name);
			columnNames.add(QSPRConstants.DM + name + ":" + QSPRConstants.CONSENSUS_STD);
		}

		model.setColumnsNames(columnNames);

		logger.info("Finished construction of the consensus model for each property");

		model.attachment.updateObject();
		model.updateDescription();
		Globals.session().saveOrUpdate(model);

		for (ModelMapping mm : model.modelMappings) 
		{
			ModelStatistics stats = (ModelStatistics) mm.statisticsRecalculated.getObject();
			stats.validationSetId="all"; // to save all sets
			stats.actualizeStatistics(mm);
			stats.recalculateStatistics(mm);
			mm.statisticsRecalculated.setObject(stats);
			mm.statisticsOriginal.updateObject();
			Globals.session().saveOrUpdate(mm);
		}
		setStatus("Finished");
	}

	private void selectBestAverage(int property, List<Model> initialChildModels, ApplicabilityDomain[][] ads, final ConsensusModelConfiguration consConf) throws UnitConversionException, IOException {

		logger.info("Creating empty statistics for property " + (property +1));

		List<Model> childModels = new ArrayList<Model>(initialChildModels); // this list can be changed...

		ModelMapping mm = model.modelMappings.get(property);
		ModelStatistics ms = ModelStatistics.getEmptyStatistics(mm);
		ms.setBootstrapReplicas(QSPRConstants.NO_REPLICAS); // to make calculation faster

		logger.info("Averaging for property " + (property +1));

		// Fix for implicit values
		Model oneChildModel = childModels.get(0);
		ModelStatistics oneChildStats = (ModelStatistics) oneChildModel.modelMappings.get(property).statisticsRecalculated.getObject(); // cached object
		SetStatistics oneTrainChildSetStats = oneChildStats.getSetStatisticsByBasket(ms.sets.get(0).basketId);

		if(oneTrainChildSetStats.hasVirtual())
			for(PointStatistics ps:oneTrainChildSetStats.points) 
				if(ps.virtual != null) 
					ms.sets.get(0).points.add(new PointStatistics(ps));

		if(consConf.classes ==null || consConf.classes.size() <= property) {
			if(consConf.classes == null)consConf.classes= new ArrayList<Integer>();
			consConf.classes.add(oneTrainChildSetStats.getNumberOfClasses());
		}

		logger.info("Averaging for property " + (property +1)+ " which has " + consConf.classes.get(property));

		// Fix for implicit values stops

		Map<Integer,Float> hashedValues = null;

		if(consConf.type == ConsensusType.OPTIMAL) {
			List<SetStatistics> allsets = new ArrayList<SetStatistics>();
			allsets.addAll(ms.sets);
			Collections.sort(childModels, new SortByPerformance(property));
			List<Model> newmodels = new ArrayList<Model>(); 
			List<Model> bestmodels = new ArrayList<Model>(); 
			double error = Double.MAX_VALUE;
			List<SetStatistics> sets = new ArrayList<SetStatistics>();
			SetStatistics ss= ms.sets.get(0);
			sets.add(ss);
			hashedValues = new HashMap<Integer,Float>();
			for(int i=0;i<childModels.size();i++) {
				if(newmodels.contains(childModels.get(i)))continue; // to avoid to have the same model two times
				newmodels.add(childModels.get(i));
				consensusForProperty(property, sets, mm.unit, newmodels, ads, consConf, hashedValues);
				if(i == 0)ms.actualizeStatistics(mm); // ??? only once since set is the same!
				ms.recalculateStatistics(mm);
				double performance = sets.get(0).getPerformanceError();
				if(performance < error*0.995) { // skip models with less than 1% difference and preferring models with smaller number of models
					error = performance;
					bestmodels.clear();
					bestmodels.addAll(newmodels);
					logger.info("property: " + property + " added new model with error: " + error + " total: " + bestmodels.size());
				}else {
					newmodels.remove(childModels.get(i)); // maybe another model will be better ?
					logger.info("property: " + property + " nothing changed :" + performance +" total: " + newmodels.size() );
				}
			}
			childModels.clear();
			childModels.addAll(bestmodels);
			ArrayList<Long> models = new ArrayList<Long>(childModels.size());
			for(Model m:childModels)
				models.add(m.publicId);
			consConf.addModelsForProperty(models);
			ms.sets = allsets; // some sets are sometimes deleted ....
			logger.info("property: " + property + " selected: " +  bestmodels.size() + " best models");
		} 

		consensusForProperty(property, ms.sets, mm.unit, childModels, ads, consConf,hashedValues);
		ms.containsPredictions = true;
		mm.statisticsOriginal = new Attachment<ModelStatistics>(ms, AttachmentSource.ModelStatistics);
		mm.statisticsRecalculated = new Attachment<ModelStatistics>(ms, AttachmentSource.ModelStatistics);
	}

	private void consensusForProperty(int property, final List<SetStatistics> sets, Unit unit, List<Model> childModels, final ApplicabilityDomain[][] ads, final ConsensusModelConfiguration consConf, Map<Integer,Float> hashedValues) throws UnitConversionException, IOException {

		Map<String, Map<Long,Integer>> remap = new HashMap<String, Map<Long,Integer>>();

		for (SetStatistics ss : sets) // going by individual sets
		{
			if(!ss.distancesToModel.contains(QSPRConstants.CONSENSUS_STD))
				ss.distancesToModel.add(QSPRConstants.CONSENSUS_STD);
			for (int psNum = 0; psNum < ss.points.size(); psNum++) // going by all predictions (rows)
			{
				PointStatistics ps = ss.points.get(psNum);
				PointStatistics pps = averageValues(ps.id, psNum, property, unit, ss, childModels, ads, remap, consConf, hashedValues);
				ps.predicted = pps.predicted;
				ps.distancesToModel = pps.distancesToModel;
				ps.error = pps.error;
				ps.ensemblePredictions = pps.ensemblePredictions;
			}
		}
	}

	/**
	 * 
	 * @param recordId - for which averaging is performed
	 * @param psNum -- data point for which averaging is performed
	 * @param property 
	 * @param unit 
	 * @param ss - set statistics
	 * @param childModels
	 * @param ads
	 * @param remap
	 * @param config
	 * @param hashedValues
	 * @return
	 * @throws UnitConversionException
	 * @throws IOException
	 */
	private PointStatistics averageValues(long recordId, int psNum, int property, Unit unit, SetStatistics ss, List<Model> childModels, ApplicabilityDomain[][] ads, Map<String, Map<Long,Integer>> remap, final ConsensusModelConfiguration config, Map<Integer,Float> hashedValues) throws UnitConversionException, IOException {

		float predictedValues[] = new float[childModels.size()];
		float predictedErrors[] = new float[childModels.size()];
		String error = null;

		for	(int i = 0; i < childModels.size(); i++)  // going by individual models for the same property
		{
			Model childModel = childModels.get(i);

			int hash = ("" + childModel.id + "_" + recordId + "_val").hashCode();

			if(hashedValues != null && hashedValues.containsKey(hash)){
				predictedValues[i] = hashedValues.get(hash);
				predictedErrors[i] = hashedValues.get(("" + childModel.id + "_" + recordId + "_err").hashCode());
				continue;
			}

			ModelStatistics childStats = (ModelStatistics) childModel.modelMappings.get(property).statisticsRecalculated.getObject(); // cached object
			SetStatistics childSetStats; 
			if (ss.basketId != null && !ss.setId.equals(QSPRConstants.EXCLUDED))
				childSetStats = childStats.getSetStatisticsByBasket(ss.basketId);
			else
				childSetStats = childStats.getSetStatistics(ss.setId);
			if (childSetStats == null)
			{
				Basket unknownBasket = Basket.getBasket(Globals.userSession(), ss.basketId);
				throw new UserFriendlyException("Model " + childModel.name + " is incompatible with the consensus, it has no data for set for a basket " + (unknownBasket == null ? ss.basketId : unknownBasket.name));
			}

			// It is also defined by units of ConsensusConfiguration
			boolean classificationProperty = childSetStats.classificationSummary != null;

			//if(classificationProperty && childSetStats.classificationSummary.getMaxClass() > 2)
			//	throw new UserFriendlyException("Consensus model currently cannot be used for properties > 2 classes. In the selected models you have property with "+childSetStats.classificationSummary.getMaxClass()+" classes");

			PointStatistics childPS = childSetStats.points.get(psNum);
			if (childPS.id != recordId){
				String key = ss.basketId+"_"+i+"_"+ss.setId.equals(QSPRConstants.EXCLUDED);
				if(!remap.containsKey(key))
					remap.put(key, createMap(childSetStats.points)); // create and cache remapping for each model
				Integer newId = remap.get(key).get(recordId);
				if(newId != null)childPS = childSetStats.points.get(newId);
				if (childPS == null)
					logger.info("Can't resolve childPS.id=" + recordId);
				else
					if (childPS.id != recordId)							
						throw new UserFriendlyException("Model <" + childModel.name + "> is incompatible with concensus (Set "+ss.setId+": "+ss.points.size()+" in the consensus model and "+childSetStats.points.size()+" in the sub-model). The sub-model must be recalculated.");
			}

			if(childPS == null || childPS.error != null){
				predictedValues[i] = predictedErrors[i] = Float.NaN;  // local values 
				error = "Error in a sub-model " + childModel.name + ": " + 
						(childPS == null?"This record was not present in the set: your training set was incompatible with the current version of the basket": childPS.error);
			}
			else{
				predictedValues[i] = (float) UnitConversion.convert(childPS.predicted, childModel.modelMappings.get(property).unit, unit, Repository.record.getWeight(recordId));

				switch(config.type){
				case BEST_MODEL:

				case WEIGHTED_AVERAGE:
					predictedErrors[i] = (float) ads[i][property].getPredictedError(childPS.distancesToModel.get(0)); //RMSE is squared for correct weighting. TODO: use the same units !!!!
					if(classificationProperty)predictedErrors[i] = 1 - predictedErrors[i]; // converting accuracy to error						
					break;
				case AVERAGE:
				case OPTIMAL:
					break;
				case RMSE_WEIGHTED:
					if(config.weights.get(i)[property] == null) // first set -- determining weighting
						config.weights.get(i)[property] = classificationProperty? 
								1 - (float)childSetStats.classificationSummary.auc.getValue(): //by AUC: any better ideas?
									(float) childSetStats.rmse ;
					predictedErrors[i] = config.weights.get(i)[property];
					break;
				}
			}

			if(hashedValues != null) {
				hashedValues.put(hash,predictedValues[i]);
				hashedValues.put(("" + childModel.id + "_" + recordId + "_err").hashCode(),predictedErrors[i]);
			}
		}

		PointStatistics ps = new PointStatistics();
		ps.id = recordId;

		AveragingAggregator agregator = new AveragingAggregator(1, config.type, config.allowErrors, true); // just for one point and for one property
		agregator.stdDmName = QSPRConstants.CONSENSUS_STD;

		if(config.classes != null && config.classes.get(property)>2)  // fooling the Aggregator
			AveragingAggregator.substituteWithMaxClass(predictedValues, predictedErrors, config.classes.get(property));

		agregator.addPredictions(0, predictedValues, predictedErrors);
		DataTable tab = agregator.getAggregatedResult(false);

		if(!tab.getRow(0).isError()){
			ps.predicted = (Double)tab.getValue(0, 0);

			ps.distancesToModel.add((Double)tab.getValue(0, 1));
			ps.error = null; // overrides error in the individual models
			ps.ensemblePredictions = (float[])tab.getValue(QSPRConstants.INDIVIDUAL_PREDICTIONS);
		}else{
			ps.distancesToModel.add(0.); // required to avoid index of out range exception
			ps.error = ps.error == null ? tab.getRow(0).detailedStatus : error; // overrides error in the individual models
		}

		return ps;
	}


	private boolean compareConditions(List<ShortCondition> first, List<ShortCondition> second) {
		if(second == null && first == null) return true;

		if(second == null && first != null) return false;
		if(second != null && first == null) return false;

		if(first.size() != second.size())return false;

		for(int i =0; i<first.size(); i++) {
			ShortCondition one = first.get(i);
			boolean found = false;
			for(int j =0; j<second.size(); j++) {
				ShortCondition two = second.get(j);
				if(one.id == two.id && one.unitId == two.unitId)found = true;
			}

			if(!found) return false;
		}

		return true;
	}

	/**
	 * Required to account bugs with change of the order of entries in baskets
	 * @param points
	 * @return
	 */

	private Map<Long, Integer> createMap(List<PointStatistics> points) {
		Map <Long, Integer> map = new HashMap<Long,Integer>();

		for( int psNum = 0; psNum < points.size(); psNum++){
			map.put(points.get(psNum).id, psNum);
		}
		return map;
	}

	private void mergeSets(Model m)
	{
		logger.info("Merging sets...");
		m.trainingSet = Basket.getById(m.trainingSet.id);
		Globals.restartMainTransaction(true);

		m.trainingSet = (Basket) Globals.session().merge(m.trainingSet);
		List<Basket> valSets = new ArrayList<Basket>();
		for (Basket b : m.getValidationSets())
			valSets.add((Basket) Globals.session().merge(b));

		m.cleanValidationSets();
		for (Basket basket : valSets)
			m.addValidationSet(basket);
	}
}

class SortByPerformance implements Comparator<Model> 
{ 
	int property;

	SortByPerformance(int property){
		this.property = property;
	}

	// Used for sorting in ascending order of 
	// roll number 
	public int compare(Model a, Model b) 
	{ 
		ModelStatistics stata = (ModelStatistics) a.modelMappings.get(property).statisticsRecalculated.getObject(); // cached object
		ModelStatistics statb = (ModelStatistics) a.modelMappings.get(property).statisticsRecalculated.getObject(); // cached object
		return Double.compare(stata.sets.get(0).getPerformanceError(),statb.sets.get(0).getPerformanceError());
	} 
} 