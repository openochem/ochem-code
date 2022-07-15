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

import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;

import qspr.Globals;
import qspr.dao.Repository;
import qspr.entities.Basket;
import qspr.entities.Model;
import qspr.metaserver.configurations.ConsensusModelConfiguration;
import qspr.modelling.applier.ModelApplier;
import qspr.workflow.utils.QSPRConstants;

public class ConsensusApplierTest
{
	@Rule
	public LoginRule rule = new LoginRule("CAT", true);
	
	@Test
	public void applyConsensusTest() throws Exception {
		
		Basket trainingSet = OCHEMTestHelpers.generateRandomBasketNumeric("CAT-Train", 100);
		Model model1 = OCHEMTestHelpers.trainAModel("CAT", trainingSet, null, "classpath:tests/ann-cv.ochem.xml", true, true);
		Model model2 = OCHEMTestHelpers.trainAModel("CAT", trainingSet, null, "classpath:tests/ann-cv.ochem.xml", true, true);
		
		
		Model consensus = OCHEMTestHelpers.trainAModel("CAT", trainingSet, null, "classpath:tests/ann-cv.ochem.xml", true, true);
		consensus.template = Repository.modelTemplate.getByName(QSPRConstants.CONSENSUS);
		ConsensusModelConfiguration consConf = new ConsensusModelConfiguration();
		consConf.addModel(model1.publicId);
		consConf.addModel(model2.publicId);
		consensus.attachment.getObject().configuration = consConf;
		consensus.attachment.updateObject();
		Globals.session().saveOrUpdate(consensus.attachment);
		Globals.session().saveOrUpdate(consensus);
		Globals.restartAllTransactions(true);
		
		ModelApplier applier = new ModelApplier();
		applier.addModel(consensus);
		
		applier.compoundsProvider.basket.addMolecule("CCCCC");
		applier.startAndWait();
		
		Assert.assertFalse(applier.isError());
	}
}
