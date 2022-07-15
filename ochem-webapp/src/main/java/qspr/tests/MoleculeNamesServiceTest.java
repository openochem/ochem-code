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

import java.util.List;

import org.junit.Assert;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TestRule;
import org.junit.rules.Timeout;

import qspr.Globals;
import qspr.entities.AbstractMoleculeName;
import qspr.entities.ExperimentalProperty;
import qspr.entities.MoleculeName;
import qspr.entities.Predicate;
import qspr.entities.Property;
import qspr.entities.ValidatedFact;
import qspr.util.MoleculePeer;

import com.eadmet.business.MoleculeNamesService;
import com.eadmet.business.PaginationFilter;

@PeriodicTest(type = "general", schedule = ScheduleType.DAILY)
public class MoleculeNamesServiceTest 
{
	@Rule
	public TestRule rule = new LoginRule("MNSTest", false);

	@Rule
	public Timeout globalTimeout = new Timeout(180000);

	@Test
	@NotifyOnFailure(minConsequentFailures = 3)
	public void moleculeNameActionsTest() throws Exception
	{
		Globals.startAllTransactions();

		ExperimentalProperty ep = new ExperimentalProperty();
		ep.predicate = Predicate.get("=");
		ep.property = OCHEMTestHelpers.generateRandomPropertyCondition(Property.TYPE_NUMERIC, false);
		ep.value = 0.5;
		ep.article = OCHEMTestHelpers.generateRandomArticle();
		ep.molecule = MoleculePeer.getMolecule("O=C(Oc1ccccc1C(=O)O)C");
		ep.owner = ep.introducer = Globals.userSession().user;
		ep.unit = ep.property.defaultUnit;
		ep.resolveConflictsAndSave();

		MoleculeNamesService service = new MoleculeNamesService();
		service.setNamesForEP(ep, new String[]{"propane", "benzene", "aspirine"});
		Assert.assertEquals(3, ep.moleculenames.size());

		Globals.session().saveOrUpdate(ep);
		Globals.restartMainTransaction(true);
		ep = (ExperimentalProperty)Globals.session().get(ExperimentalProperty.class, ep.id);

		try
		{
			List<AbstractMoleculeName> names;

			service.checkNames(ep);

			Globals.session().saveOrUpdate(ep);
			Globals.restartMainTransaction(true);
			ep = (ExperimentalProperty)Globals.session().get(ExperimentalProperty.class, ep.id);

			names = service.listMoleculeNames(ep, null); //Quick check with existing molecule via VIEW
			for (AbstractMoleculeName name : names)
				if ("aspirine".equals(name.getName().toLowerCase()))
					Assert.assertEquals(1, name.getValidation());
				else
					Assert.assertEquals(4, name.getValidation());

			names = service.listMoleculeNames(ep, MoleculePeer.getMolecule("C1=CC=CC=C1")); //Slow check with candidate molecule via one by one validation

			for (AbstractMoleculeName name : names)
				if ("benzene".equals(name.getName().toLowerCase()))
					Assert.assertEquals(1, name.getValidation());
				else
					Assert.assertEquals(4, name.getValidation());

			MoleculeName mn = null;
			for(int i=0;i<ep.moleculenames.size();i++){			
				mn = ep.moleculenames.get(i);
				if(mn.name.toLowerCase().equals("benzene"))break;
			}
			service.validate(mn, ep.molecule);
			Assert.assertEquals(1, mn.getValidation());
			Assert.assertEquals(Globals.userSession().user.id.longValue(), Integer.valueOf(mn.validatedFacts.get(0).sourceid).longValue());

			Globals.restartMainTransaction(true);
			ep = (ExperimentalProperty)Globals.session().get(ExperimentalProperty.class, ep.id);
			mn = ep.moleculenames.get(0);

			service.invalidate(mn);
			Assert.assertEquals(0, mn.getValidation());


			List<MoleculeName> synonyms = service.listSynonyms(MoleculePeer.getMolecule("C1=CC=CC=C1"), null, new PaginationFilter(1, 10));

			for (MoleculeName molname : ep.moleculenames)
				Assert.assertTrue(synonyms.contains(molname)); //Our record names should appear as synonyms in this query

		} finally
		{			
			for (MoleculeName mn : ep.moleculenames)
			{
				for (ValidatedFact vf : mn.validatedFacts)
					Globals.session().delete(vf);
				mn.validatedFacts.clear();
			}
			Globals.commitAllTransactions();
		}
	}
}
