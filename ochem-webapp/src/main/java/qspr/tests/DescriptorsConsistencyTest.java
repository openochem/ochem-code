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

import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import qspr.Environment;
import qspr.Globals;
import qspr.OCHEMConfiguration;
import qspr.ThreadScope;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.Molecule;
import qspr.metaserver.configurations.DescriptorsCDDDConfiguration;
import qspr.metaserver.configurations.DescriptorsCDK2Configuration;
import qspr.metaserver.configurations.BalloonConfiguration;
import qspr.metaserver.configurations.DescriptorsAbstractConfiguration;
import qspr.metaserver.configurations.DescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsEPAConfiguration;
import qspr.metaserver.configurations.DescriptorsFragmentorConfiguration;
import qspr.metaserver.configurations.DescriptorsGSFragConfiguration;
import qspr.metaserver.configurations.DescriptorsInductiveDescriptorsConfiguration;
import qspr.metaserver.configurations.DescriptorsKrakenXConfiguration;
import qspr.metaserver.configurations.DescriptorsMAP4Configuration;
import qspr.metaserver.configurations.DescriptorsMOLD2Configuration;
import qspr.metaserver.configurations.DescriptorsMersyConfiguration;
import qspr.metaserver.configurations.DescriptorsMOPAC2016Configuration;
import qspr.metaserver.configurations.DescriptorsMOPACConfiguration;
import qspr.metaserver.configurations.DescriptorsMORDREDConfiguration;
import qspr.metaserver.configurations.DescriptorsMeraConfiguration;
import qspr.metaserver.configurations.DescriptorsMolPrintConfiguration;
import qspr.metaserver.configurations.DescriptorsOEstateConfiguration;
import qspr.metaserver.configurations.DescriptorsPaDEL2Configuration;
import qspr.metaserver.configurations.DescriptorsPyDescriptorConfiguration;
import qspr.metaserver.configurations.DescriptorsQNPRConfiguration;
import qspr.metaserver.configurations.DescriptorsRDKITConfiguration;
import qspr.metaserver.configurations.DescriptorsSIGMAConfiguration;
import qspr.metaserver.configurations.DescriptorsSIRMSConfiguration;
import qspr.metaserver.configurations.DescriptorsSpectrophoresConfiguration;
import qspr.metaserver.configurations.StandartizationOptions;
import qspr.metaserver.configurations.DescriptorsStructuralAlertsConfiguration;
import qspr.metaserver.protocol.Task.TaskPriority;
import qspr.modelling.ModelProcessor;
import qspr.util.WrapperThread;
import qspr.workflow.datatypes.AbstractDataRow;
import qspr.workflow.datatypes.DataTable;
import qspr.workflow.datatypes.WorkflowNodeData;

import com.eadmet.business.DescriptorsCalculatorProcessor;
import com.eadmet.descriptorcache.DescriptorConfigEntry;
import com.eadmet.descriptorcache.DescriptorsCache;
import com.eadmet.utils.NumericalValueStandardizer;

/**
 * Check that descriptors calculations do not change.
 * This test used previously cached descriptors and compares them with the newly calculated predictions.
 * 
 * @author itetko
 *
 */
@RunWith(value = Parameterized.class)
@PeriodicTest(schedule = ScheduleType.DAILY, type = "general")
public class DescriptorsConsistencyTest 
{
	private static final String PREFIX_CONSTEST = "TCons";
	private static final Logger logger = LogManager.getLogger(DescriptorsConsistencyTest.class);

	/**
	 * The tolerated percentage of inconsistencies 
	 * (e.g., 10 discrepancies out of 100 predictions not larger than 10% of values are ok, may be caused by different representations)
	 */
	private static final double INCONSISTENCIES_TOLERANCE = 0.01; // ratio
	protected static final int MAXTRIALS = 2;

	private DescriptorsAbstractConfiguration configuration;
	private double difference = 0;
	static Basket basket;
	static WorkflowNodeData workflow;

	public DescriptorsConsistencyTest(DescriptorsAbstractConfiguration config, String descriptorsName)
	{
		this.configuration = config;
		logger.info("model consistency test for " + descriptorsName + " (" + config + ") constructed ");
	}

	// take the second parameter as name, in this case the model name
	@Parameters(name="{1}")
	public static Collection<Object[]> descriptors() throws Exception
	{
		Globals.startAllTransactions();

		List<DescriptorsAbstractConfiguration> descriptors = new ArrayList<DescriptorsAbstractConfiguration>();

		boolean lighweight = OCHEMConfiguration.verboseMode <= 1;

		if(OCHEMConfiguration.verboseMode >= 0) {
			descriptors.add(new DescriptorsMOLD2Configuration());
			descriptors.add(new DescriptorsMORDREDConfiguration());
			descriptors.add(new DescriptorsCDK2Configuration());
			descriptors.add(new DescriptorsOEstateConfiguration());
			descriptors.add(new DescriptorsQNPRConfiguration());
			descriptors.add(new DescriptorsInductiveDescriptorsConfiguration());
			descriptors.add(new DescriptorsSIRMSConfiguration(1,4));
			descriptors.add(new DescriptorsRDKITConfiguration());
			descriptors.add(new DescriptorsPaDEL2Configuration());
			DescriptorsStructuralAlertsConfiguration alerts =  Repository.various.getEFGConfiguration();
			if(alerts.alertPatterns != null && alerts.alertPatterns.size() > 0)
				descriptors.add(alerts);
			try {
				descriptors.add((DescriptorsAbstractConfiguration) Class.forName("qspr.metaserver.configurations.DescriptorsEPAConfiguration").newInstance());
			} catch (ClassNotFoundException e) {
				logger.warn("Skipping tests for ChemaxonScaffold. Class was not found. Chemxon plugins not loaded?");
			}
			descriptors.add(new DescriptorsSpectrophoresConfiguration());

			descriptors.add(new DescriptorsGSFragConfiguration());
			descriptors.add(new DescriptorsMersyConfiguration());
			descriptors.add(new DescriptorsMeraConfiguration());
			descriptors.add(new DescriptorsFragmentorConfiguration());
			descriptors.add(new DescriptorsMolPrintConfiguration());

			descriptors.add(new DescriptorsCDDDConfiguration());
			descriptors.add(new DescriptorsMAP4Configuration());
			descriptors.add(new DescriptorsPaDEL2Configuration());
			descriptors.add(new DescriptorsEPAConfiguration());

			descriptors.add(new DescriptorsPyDescriptorConfiguration());
			descriptors.add(new DescriptorsMOPACConfiguration());
			descriptors.add(new DescriptorsMOPAC2016Configuration());
			descriptors.add(new DescriptorsKrakenXConfiguration());
			descriptors.add(new DescriptorsSIGMAConfiguration());

			try {
				descriptors.add((DescriptorsAbstractConfiguration) Class.forName("qspr.metaserver.configurations.DescriptorsChemaxonScaffoldConfiguration").newInstance());
			} catch (ClassNotFoundException e) {
				logger.warn("Skipping tests for ChemaxonScaffold. Class was not found. Chemxon plugins not loaded?");
			}
		}	

		Object[][] data = new Object[descriptors.size()][2];
		for (int i = 0; i < descriptors.size(); i++)
		{
			String name =  descriptors.get(i).getDefaultTypeName();
			try{ // we try to set all descriptors ON for all configurations, in case if they support it
				DescriptorsAbstractConfiguration conf = descriptors.get(i);
				conf.setAllOn();
				String add = ""+ conf;
				if(!name.equals(add)) {
					name = ModelConsistencyTest.normalizeOnlyLettersNumbers(name + "_" + add);
					if(name.length()>ModelConsistencyTest.LENGTH) name = name.substring(0, ModelConsistencyTest.LENGTH);
					name = name.trim();
				}
				data[i] = new Object[]{conf,name};
			}catch(Exception e){
				data[i] = new Object[]{descriptors.get(i),name};
			}
		}

		basket = getBasket(lighweight ?  "training-set.csv" : "models-set.csv" );

		DataTable mols = ModelProcessor.basketToSDFTable(basket);
		workflow = new WorkflowNodeData();
		workflow.ports.add(mols);

		Globals.rollbackAllTransactions();

		return Arrays.asList(data);
	}

	/**
	 * Test all the featured models (classification and regression models) defined as parameters
	 * @throws Exception
	 */
	@Test(timeout = 7200000) // 2 hours per descriptors set
	public void descriptorsConsistencyTest() throws Exception
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

						try {
							PeriodicTests.checkIfTaskSupported(configuration.getDefaultTypeName());
						}catch(Exception e) {
							trial = MAXTRIALS;
							throw e;
						}

						DescriptorsCache cache = new DescriptorsCache(false,false);
						DescriptorConfigEntry entry = new DescriptorConfigEntry(configuration, configuration.getDefaultTypeName());
						DataTable oldPredictions = cache.getDescriptorsFromCache(entry, workflow);

						DescriptorsCalculatorProcessor processor = new DescriptorsCalculatorProcessor();
						processor.basket = basket;

						processor.defaultTaskPriority = TaskPriority.HIGH + 1;

						processor.structureStandartisation = new StandartizationOptions(true);

						configuration.moleculeTimeout = 3; //three minutes per molecule maximum! some calculations may fail because of timeout

						if(configuration.requires3D())
							processor.structureOptimisation = new BalloonConfiguration();
						processor.descConfig = new DescriptorsConfiguration();
						processor.descConfig.addDescriptorType(configuration.getDefaultTypeName(), configuration);
						processor.descConfig.forceUpdateDescriptorCache = true; // of course we recalculate all descriptors

						Globals.rollbackAllTransactions();

						processor.start();

						while (!processor.isReady())
							Thread.sleep(1000);

						DataTable newPredictions = processor.dtDescriptors;			

						if(oldPredictions.columns == null || oldPredictions.columns.size() == 0)
							oldPredictions = cache.getDescriptorsFromCache(entry, workflow); // try again, maybe these values were already automatically cached!

						if(oldPredictions.columns == null || oldPredictions.columns.size() == 0 || oldPredictions.getRowsNoErrorsSize() < oldPredictions.getRowsSize()/2){ // still not, let us save them
							if(entry.user != null)throw new AssertionError("User is not null: " + entry.user);
							cache.saveNewValues(entry, workflow, newPredictions);
							throw new  Exception("No descriptors were found - recaching");
						}
						else{
							messages = comparePredictions(newPredictions, oldPredictions, configuration.getDefaultTypeName().contains("MOPAC"));

							if (!messages.isEmpty()) {
								String message = configuration.getDefaultTypeName() + " had " + messages.size() + " inconsistent descriptors out of " + oldPredictions.getRowsSize() + 
										" \nExample:" + messages.toArray()[messages.size()-1];
								if (messages.size() > oldPredictions.getRowsSize() * INCONSISTENCIES_TOLERANCE )
									throw new Exception(message);
								else
									logger.info(message);
							}
							trial = MAXTRIALS;
						}
					}catch(Exception e){
						Thread.sleep(10000); // to allow predictions to cache
						if(trial >= (MAXTRIALS-1)) 
							throw new AssertionError(e.getMessage() == null?e:e.getMessage());
						else
							logger.info((e.getMessage() == null?e:e.getMessage()) + " - starting again!");
					}
			}
		};

		wrapper.run();

		if (wrapper.exception != null)
			throw wrapper.exception;
	}


	private Collection<String> comparePredictions(DataTable newPredictions, DataTable cachedPredictions, boolean ignoreFailure) throws Exception
	{
		Map<Double,String> messages = new TreeMap<Double,String>();

		if(newPredictions.getRowsSize() != cachedPredictions.getRowsSize())
			throw new Exception("Inconsistent number of records");

		double all=0;
		difference = 0;

		Map<Integer,Integer> newToCached = new HashMap<Integer,Integer>();

		for(int i=0; i<newPredictions.columns.size(); i++)
			newToCached.put(i, cachedPredictions.columns.indexOf(newPredictions.columns.get(i)));

		int size = newPredictions.getRowsSize();

		for (int i = 0; i < size; i++)
		{
			if(cachedPredictions.getRow(i).isError() || DescriptorsCache.SENT_FOR_CALCULATIONS.equals(cachedPredictions.getRow(i).status))continue; // old errors skipping this molecule

			AbstractDataRow row = newPredictions.getRow(i);

			if (row.isError()){
				String status = row.detailedStatus == null ? "" : row.detailedStatus.toLowerCase();
				if(!status.contains("timeout") && !ignoreFailure ) {
					messages.put(Double.POSITIVE_INFINITY, " Failed record n = " + i + " " + row.detailedStatus);
				}
				else{
					status = "Molecule mapping2: " + basket.entries.get(i).ep.molecule.mapping2.id + " n = " + i + " failed with timeout";
					messages.put(1000., status);
					logger.error("DESCRIPTORS: " + status);
					continue;
				}
			}

			boolean first = true;

			for(int j=0; j < newPredictions.getColumnsSize();j++){ // only prediction results columns
				int oldi = newToCached.get(j);
				double newRes = NumericalValueStandardizer.getSignificantDigitsDouble(Double.valueOf(""+ newPredictions.getValue(i,j)), 3);
				double oldRes = oldi >= 0? NumericalValueStandardizer.getSignificantDigitsDouble(Double.valueOf(""+cachedPredictions.getValue(i,oldi)),3):0; // some columns are always 0 and may not be saved

				double diff  = 0, average = 0;

				if(Double.isFinite(newRes) || Double.isFinite(newRes)) { // if we have both times NaN - this is OK.
					diff = Math.abs(newRes-oldRes);
					average = Math.abs(newRes+oldRes)/2.;
				}

				difference += diff;all += average;
				if(diff == 0 || diff < average*INCONSISTENCIES_TOLERANCE)continue;

				if(first){   // only one message per molecule
					String str = " new prediction " + newRes + " != " + oldRes + " for basket entry: " + i + " descriptor: " + newPredictions.getColumn(j);
					messages.put(Math.abs(diff), str); 
					logger.info("DESCRIPTORS: " + configuration + str);
					first = false;
				}
			}
		}

		difference =  Math.round(10000*difference/all)/100;
		if(difference > 0 && messages.size() > 0)logger.info("Difference: " + difference + "% for " + configuration + " warnings: " + messages.size());
		return messages.values();
	}


	static public Basket getBasket(String fileName) throws Exception{
		File f = new File(getTestsFolder()+"/"+fileName);

		logger.info("Using file: "+f.getCanonicalPath());

		MoleculesWithValues mols = MoleculesWithValues.loadFromFile(new FileInputStream(f));

		Basket basket = new Basket();
		for (Molecule mol : mols.molecules)
			basket.addMolecule(mol);

		return basket;
	}


	static protected String getTestsFolder(){
		if (Environment.rootDir == null || Environment.rootDir.equals(""))
			Environment.rootDir = new File(ModelCreationTest.class.getClassLoader().getResource("views.properties").getPath()).getParentFile().getAbsolutePath() + "/../../src/main/webapp";
		return Environment.getWebInfDir() + "/tests/models";
	}

}