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

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.RuleChain;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import qspr.Globals;
import qspr.ThreadScope;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.modelling.applier.ModelApplier;
import qspr.modelling.applier.ModelApplierTaskProcessor;

/**
 * 
 * Run all the published models and check whether they are still valid.
 * 
 * This test only checks the validity (whether the models run at all) but not consistency.
 * It might be a good idea to use the consistency test for all the published models instead.
 *
 * @author midnighter
 *
 */
@PeriodicTest(schedule = ScheduleType.AT_FIVE_AM, type = "general")
@SuppressWarnings("unchecked")
public class PublishedModelsTest
{

	//////// Currently disabled 

	@Rule													
	public TestRule chain = RuleChain.outerRule(new Timeout(3600000))
	.around(new LoginRule("PublModels", false));

	@Test
	@NotifyOnFailure(minConsequentFailures = 2)
	public void runAllPublished() throws Exception 
	{
		for(int i=0; i<2; i++){

			if(i<3)continue; //  disabled sine we use much more powerful ones //TODO remove if required

			ModelApplier applier = new ModelApplier();

			logger.info("ThreadScope.get().userSession: " + ThreadScope.get().userSession.id);
			ThreadScope.get().userSession.disableQuota = true;

			Globals.startAllTransactions();
			List<Model> models = Globals.session().createQuery("from Model where published=1 and taskId is null and approved=1 " 
					+ "and configurationXml "+(i==1?"not":"")+" like '%cons%' and length(configurationXml) > 0 ").list();

			logger.info("Testing models: "+models.size());

			Collections.sort(models, new CompareModels());  

			Set<Long> ignored = new HashSet<Long>(Arrays.asList(484l));

			for (Model model : models)
			{
				if(ignored.contains(model.publicId))continue;
				logger.info("Testing "+model.name + " size=" + model.size);
				applier.addModel(model);
			}

			applier.compoundsProvider.basket = new Basket();
			applier.compoundsProvider.basket.addMolecule("C=C=CCC=CCCCCCC=C=C=C=C=C=C=CC(=C)CC");
			applier.useCache = false;

			applier.start();
			Globals.commitAllTransactions();

			logger.info("Starting waiting...");
			applier.awaitTasksFinished();
			logger.info("Finished waiting...");

			int failedModels = 0;
			for (ModelApplierTaskProcessor mTask : applier.modelTasks)
				if (mTask.isError())
				{
					logger.error(mTask.model.name + "  mid=" + mTask.model.publicId+" has Failed with error " + mTask.getStatus() + ". ");
					failedModels++;
				}

			logger.info(String.format("%d out of %d models have Failed", failedModels, models.size()));

			Assert.assertEquals(0, failedModels);
		}
	}


	private static final Logger logger = LogManager.getLogger(PublishedModelsTest.class);
}

class CompareModels implements Comparator<Model>{

	@Override
	public int compare(Model o1, Model o2) {
		long v = o2.size-o1.size;
		return v > Integer.MAX_VALUE? Integer.MAX_VALUE : v < Integer.MIN_VALUE? Integer.MIN_VALUE: (int) v;
	}

}
