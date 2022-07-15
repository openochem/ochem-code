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

package qspr.tests;

import java.io.StringReader;
import java.text.Normalizer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.Unmarshaller;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.configuration.JAXBContextFactory;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.BasketEntry;
import qspr.entities.Model;
import qspr.entities.ModelAttachment;
import qspr.metaserver.configurations.DescriptorType;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.modelling.ModelApplierCache;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.configurations.CDSConfiguration;
import qspr.util.WrapperThread;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.utils.QSPRConstants;

import com.eadmet.utils.NumericalValueStandardizer;

/**
 * Check that model predictions do not change.
 * This test used previously cached predictions and compares them with the newly calculated predictions.
 * The test is running for the published models that are of especial importance for us.
 * 
 * Some percentage of inconsistencies is tolerated
 * 
 * @author midnighter/itetko
 *
 */
@RunWith(value = Parameterized.class)
@PeriodicTest(schedule = ScheduleType.DAILY, type = "general")
public class ModelConsistencyTest 
{
	private static final String PREFIX_CONSTEST = "MCons";
	private static final Logger logger = LogManager.getLogger(ModelConsistencyTest.class);

	/**
	 * The tolerated percentage of inconsistencies 
	 * (e.g., 10 discrepancies out of 100 predictions not larger than 10% of values are ok, may be caused by different representations)
	 */
	private static final double INCONSISTENCIES_TOLERANCE = 0.01; // fraction
	protected static final int MAXTRIALS = 2;

	public static int LENGTH = 60;
	private Long pubId;
	private double difference = 0;
	static Basket basket;

	public ModelConsistencyTest(Long modelId, String modelName)
	{
		this.pubId = modelId;
		logger.info("model consistency test for " + modelName + " (id: " + modelId + ") constructed ");
	}

	// take the second parameter as name, in this case the model name
	@Parameters(name="{1}")
	public static Collection<Object[]> models() throws Exception
	{

		Globals.startAllTransactions();

		boolean cleanCache = false;

		List<Model> all = Repository.model.getPublishedApprovedModels();
		if(OCHEMConfiguration.verboseMode == -1) {
			int ids[] = {9, 613}; 
			all.clear(); for(int l: ids)all.add(Repository.model.getByPublicId(l));
		}

		boolean lighweight = OCHEMConfiguration.verboseMode <= 1;

		List<Model> models = new ArrayList<Model>();

		JAXBContext context = JAXBContextFactory.get("qspr.modelling.configurations:qspr.metaserver.configurations:qspr.entities");
		Unmarshaller unmarshaller = context.createUnmarshaller();
		Set<DescriptorType> types = new HashSet<DescriptorType>();

		for(Model model:all){
			if(QSPRConstants.UPLOADED_MODEL.equals(model.template.name))continue; // no test for uploaded models
			if(model.name.contains("McReynolds"))continue; 
			if(lighweight && !analyse(model.publicId,unmarshaller,types))continue;
			models.add(model);
			if(cleanCache) model.clearCache();  // to delete previous cache
		}

		basket=DescriptorsConsistencyTest.getBasket(lighweight ?  "training-set.csv" : "models-set.csv");

		Globals.commitAllTransactions();

		int trials = OCHEMConfiguration.verboseMode == -1? 10 :1;

		Object[][] data = new Object[models.size()*trials][2];
		for (int i = 0; i < models.size(); i++)
		{
			for(int j =0; j<trials;j++) {
				String name =  models.get(i).name;
				name = Normalizer.normalize(name, Normalizer.Form.NFD);
				name = normalizeOnlyLettersNumbers(name);
				name = j == 0? name:name + "_"+j;
				if(name.length()>LENGTH)name=name.substring(0, LENGTH);
				name = name +  "_" + models.get(i).publicId;
				data[i*trials +j] = new Object[]{models.get(i).publicId,name};
			}
		}

		logger.info(models.size() + " will be sent for consistency check, ignored " + (all.size() - models.size()));

		return Arrays.asList(data);
	}

	static public boolean isBlank(String value) {
		return (value == null || value.equals("") || value.equals("null") || value.trim().equals(""));
	}

	static public String normalizeOnlyLettersNumbers(String str) {
		String s;
		if (!isBlank(str)) {
			s=str.replaceAll("[^\\p{L}\\p{Nd}-]+", " ");
		} else {
			s = "";
		}
		s.replaceAll("\\s+", " ");
		return s;
	}

	/**
	 * Test all the featured models (classification and regression models) defined as parameters
	 * @throws Exception
	 */
	@Test(timeout = 1800000) // 30 min per model
	public void modelConsistencyTest() throws Exception
	{
		WrapperThread wrapper = new DataDrivenTestWrapper(PREFIX_CONSTEST) 
		{
			@Override
			public void wrappedTest() throws Exception 
			{
				Collection<String> messages = null;
				for(int trial = 0; trial < MAXTRIALS; trial++)
					try{
						Globals.restartAllTransactions(true);
						ThreadScope.get().userSession.disableQuota = true;

						Model model = Repository.model.getByPublicId(pubId);

						logger.info("Starting model " + model.name + " id: " + model.publicId);

						ModelApplier applier = new ModelApplier();
						applier.addModel(model);
						applier.useCache = false;
						applier.defaultTaskPriority = TaskPriority.HIGH + 1;

						DataTable oldPredictions = getMoleculesAndCachedPredictions(model); //getting them from cache using another applier
						int oldrecords = oldPredictions.getRowsNoErrorsSize();

						if (oldPredictions.getRowsNoErrorsSize() > oldPredictions.getRowsSize()/2)
							logger.info("Using " + oldPredictions.getRowsNoErrorsSize() + " cached predictions to validate model consistency");
						else{
							logger.info("Re-caching predictions for new analyses");
							applier.useCache = true;
						}

						applier.compoundsProvider.basket = basket;

						logger.info("Using " + applier.compoundsProvider.basket.getRowsSize() + " out of " 
								+ oldPredictions.getRowsSize() + " cached predictions to validate model consistency");

						applier.start();

						while (!applier.isReady())
						{
							Thread.sleep(500);
							applier.update();
						}

						if (applier.isError())
							throw new Exception("Applier failed: " + applier.getErrorMessage());

						if(applier.modelTasks.get(0).isError())
							throw new AssertionError("model: "+ pubId +" " + model.name + "\tcompletely failed");

						DataTable newPredictions = applier.modelTasks.get(0).wndResult.ports.get(0);

						if(oldrecords <  oldPredictions.getRowsSize()/2) 
							throw new Exception("Old predictions were mostly errors." + newPredictions.getRowsNoErrorsSize() + " new predictions will be cached out of : " + newPredictions.getRowsSize()); // no check at this moment of time


						messages = comparePredictions(newPredictions, oldPredictions, applier.compoundsProvider.basket.entries,model.name);

						trial = MAXTRIALS; // nothing serious! we can quit now

						if (!messages.isEmpty()) {
							String message = model.name + " ID: " + model.publicId + " had " + messages.size() + " inconsistent predictions out of " + oldPredictions.getRowsSize() + " \nExample: " + messages.toArray()[messages.size()-1];

							if (messages.size() > oldPredictions.getRowsSize() * INCONSISTENCIES_TOLERANCE )
								throw new Exception(message);
							else
								logger.info(message);
						}

					}catch(Exception e){
						Thread.sleep(10000); // to allow predictions to cache
						if(trial > 0) 
							throw new AssertionError(e.getMessage() == null?e:e.getMessage());
						else
							logger.info((e.getMessage() == null?e:e.getMessage()) + " - starting again!");
					}
			}

			/*
			private boolean oldCdkId(int pubId) {
				Integer ids[] = {29,62,63,148,303,350,360,370,380,390,400,410,420,430,440,450,460};
				for(Integer i: ids) 
					if(pubId == i) return true;
				return false;
			}
			 */		};

			 wrapper.run();

			 if (wrapper.exception != null)
				 throw wrapper.exception;
	}

	static boolean analyse(long publicId, Unmarshaller unmarshaller, Set<DescriptorType> types){

		Model model = Repository.model.getByPublicId(publicId);
		if(QSPRConstants.CONSENSUS.equals(model.template.name))return false; // no test for consensus models

		boolean analyse = false;

		try{
			String s=model.configurationXml.replace("configuration xsi", "configuration xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi");
			StringReader sr = new StringReader(s);
			ModelAttachment mod = (ModelAttachment) unmarshaller.unmarshal(sr);
			if(!(mod.configuration instanceof CDSConfiguration)) return false;
			for(DescriptorType type: ((CDSConfiguration)mod.configuration).descriptors.types){
				boolean found = false;
				for(DescriptorType typ:types)
					if(typ.similar(type))
						found = true;
				if(!found){ 
					analyse = true;
					logger.error("New descriptor type: " + type + " " + types.size() + " " +  model.name + " id: " + model.publicId);
				}
				types.add(type);
			}

			return analyse;

		}catch(Exception e){
			logger.error(e);
			logger.error(model.configurationXml);
			return true;
		}
	} 

	private Collection<String> comparePredictions(DataTable newPredictions, DataTable cachedPredictions, List<BasketEntry> entries, String name) throws Exception
	{
		Map<Double,String> messages = new TreeMap<Double,String>();

		if(newPredictions.compareColumns(cachedPredictions) != null)
			throw new Exception("Inconsistent column names: " + newPredictions.compareColumns(cachedPredictions));

		if(newPredictions.getRowsSize() != cachedPredictions.getRowsSize())
			throw new Exception("Inconsistent number of records");

		boolean found = false;
		for(int j=0; j < newPredictions.getColumnsSize();j++) // only prediction results columns
			if(newPredictions.columns.get(j).contains(QSPRConstants.PREDICTION_RESULT_COLUMN))found = true;

		if(!found && newPredictions.getRowsNoErrorsSize() > 0)throw new Exception("The predictions should contain at least one " + QSPRConstants.PREDICTION_RESULT_COLUMN + "column");

		double all=0;
		difference = 0;

		for (int i = 0; i < newPredictions.getRowsSize(); i++)
		{
			if(cachedPredictions.getRow(i).isError())continue; // old errors skipping this molecule

			if (newPredictions.getRow(i).isError()){

				double oldRes = NumericalValueStandardizer.getSignificantDigitsDouble(Double.valueOf(""+cachedPredictions.getValue(i,0)), 3);

				String str = "New prediction failed: "  + "model: " + pubId + " for record: " + i + " old value: " + oldRes +
						" for mappping2 M" + entries.get(i).ep.molecule.mapping2.id + "\t("+name+") with status = "
						+ newPredictions.getRow(i).detailedStatus + "\n\n\n"; // new error!
				messages.put(Double.MAX_VALUE, str);

				double diff = Math.abs(oldRes); if(diff == 0) diff = 1;
				difference += diff;all += diff;
				continue;
			}

			boolean first = true;

			for(int j=0; j < newPredictions.getColumnsSize();j++) // only prediction results columns
				if(newPredictions.columns.get(j).contains(QSPRConstants.PREDICTION_RESULT_COLUMN)){

					double newRes = NumericalValueStandardizer.getSignificantDigitsDouble(Double.valueOf(""+ newPredictions.getValue(i,j)), 3);
					double oldRes = NumericalValueStandardizer.getSignificantDigitsDouble(Double.valueOf(""+cachedPredictions.getValue(i,j)), 3);

					double diff = Math.abs(newRes-oldRes);
					double average = Math.abs(newRes+oldRes)/2.;
					difference += diff;all += average;
					if(diff == 0 || diff < average*INCONSISTENCIES_TOLERANCE)continue;

					if(first){  // only one message per molecule
						String str = "model: " + pubId + " new prediction for " + i + " " + newRes + " != " + oldRes + 
								" for mappping2 M" + entries.get(i).ep.molecule.mapping2.id + "\t("+name+")";
						messages.put(Math.abs(newRes - oldRes), str);
						logger.info("MODELS: "+ str);
						first = false;
					}
				}
		}

		difference  =  Math.round(10000.*difference/all)/100.;
		if(difference > 0 && messages.size() > 0)logger.info("Difference: " + difference + "% for modelid: " + pubId + " warnings: " + messages.size());
		return messages.values();
	}

	private DataTable getMoleculesAndCachedPredictions(Model model) throws Exception
	{

		ModelApplierCache cache = new ModelApplierCache(model,basket,null,null);

		DataTable cached = cache.getCachedResult("no results were found");

		Assert.assertEquals(cached.getRowsSize(), basket.entries.size());

		return cached;
	}

}